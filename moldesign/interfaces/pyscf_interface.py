import warnings
from cStringIO import StringIO
import sys

import numpy as np
from moldesign import compute

try:
    from pyscf import gto, scf, mp, mcscf, ao2mo
    import pyscf.grad
    from pyscf.dft import numint
except ImportError:
    force_remote = True
    gto = scf = mp = mcscf = ao2mo = numint = None
except OSError as exc:  # TODO: on OSX, this does ... something
    print 'WARNING: PySCF error on import - %s: %s' % (type(exc), exc)
    print 'WARNING: Using containerized version of PySCF, not the local installation'
    force_remote = True
    gto = scf = mp = mcscf = ao2mo = numint = None
else:
    force_remote = False

from moldesign.utils import if_not_none, redirect_stderr, DotDict
from moldesign import basemethods, orbitals, logs
from moldesign import units as u


# TODO: need to handle parameters for max iterations,
# level shifts, requiring convergence, restarts, initial guesses

class PySCFPotential(basemethods.QMBase):
    DEFAULT_PROPERTIES = ['potential_energy',
                          'electronic_state',
                          'mulliken']
    ALL_PROPERTIES = DEFAULT_PROPERTIES + ['eri_tensor',
                                           'forces',
                                           'nuclear_forces',
                                           'electronic_forces']

    # TODO: need to store theory name (and basis name) synonyms
    if scf is not None:  # shielding these from failed imports
        THEORIES = {'hf': scf.RHF, 'rhf': scf.RHF,
                    'uhf': scf.UHF,
                    'mcscf': mcscf.CASSCF, 'casscf': mcscf.CASSCF,
                    'casci': mcscf.CASCI,
                    'mp2': mp.MP2}
        FORCE_CALCULATORS = {'rhf': pyscf.grad.RHF, 'hf': pyscf.grad.RHF}


    NEEDS_REFERENCE = set('mcscf casscf casci mp2')
    FORCE_UNITS = u.hartree / u.bohr

    def __init__(self, **kwargs):
        super(PySCFPotential, self).__init__(**kwargs)
        self.pyscfmol = None
        self.reference = None
        self.kernel = None
        self.logs = StringIO()
        self.last_el_state = None
        self.logger = logs.Logger('PySCF interface')  # eventually, this be a child of the parent molecule

    @compute.runsremotely(remote=force_remote, is_imethod=True)
    def calculate(self, requests=None, guess=None):
        self.logger = logs.Logger('PySCF calc')  # eventually, this be a child of the parent molecule
        do_forces = 'forces' in requests
        if do_forces and self.params.theory not in self.FORCE_CALCULATORS:
            raise ValueError('Forces are only available for the following theories:'
                             ','.join(self.FORCE_CALCULATORS))
        if do_forces:
            force_calculator = self.FORCE_CALCULATORS[self.params.theory]

        if guess is not None:
            dm0 = guess.density_matrix
        elif self.last_el_state is not None:
            dm0 = self.last_el_state.density_matrix
        else:
            dm0 = None
        self.prep(force=True)  # rebuild every time
        result = {}
        if self.reference:
            kernel, failures = self._converge(self.reference, dm0=dm0)
            result['reference_energy'] = (kernel.e_tot * u.hartree).defunits()
            # TODO: These objects can't be pickled ... need to fix that.
            # result['_converged_kernel'] = kernel
            # result['_convergence_failures'] = failures
        kernel, failures = self._converge(self.kernel, dm0=dm0)
        result['potential_energy'] = (kernel.e_tot * u.hartree).defunits()
        # result['_converged_kernel'] = kernel
        # result['_convergence_failures'] = failures

        # A little bit of wfn analysis
        if self.reference is not None:
            orb_calc = self.reference
        else:
            orb_calc = self.kernel
        ao_matrices = self._get_ao_matrices(orb_calc)
        scf_matrices = self._get_scf_matrices(orb_calc, ao_matrices)
        ao_pop, atom_pop = orb_calc.mulliken_pop(verbose=-1)

        # Build the electronic state object
        basis = PySCFAOBasis(self.mol,
                             basis_fns=self._get_ao_basis_functions(),
                             h1e=ao_matrices.h1e.defunits(),
                             overlaps=ao_matrices.sao,
                             name=self.params.basis,
                             pyscfmol=self.pyscfmol)
        el_state = orbitals.ElectronicState(self.mol,
                                            self.pyscfmol.nelectron,
                                            theory=self.params.theory,
                                            aobasis=basis,
                                            ao_population=ao_pop,
                                            nuclear_repulsion=(self.pyscfmol.energy_nuc() *
                                                               u.hartree).defunits(),
                                            **scf_matrices)

        # Build and store the canonical orbitals
        cmos = []
        for coeffs, energy, occ in zip(orb_calc.mo_coeff.T,
                                       orb_calc.mo_energy * u.hartree,
                                       orb_calc.get_occ()):
            cmos.append(orbitals.Orbital(coeffs, wfn=el_state, energy=energy.defunits(), occupation=occ))
        el_state.add_orbitals(cmos, orbtype='canonical')

        # Calculate the forces
        if do_forces:
            g = force_calculator(self.kernel)
            f_e = -1.0 * g.grad_elec().reshape(self.mol.ndims) * self.FORCE_UNITS
            f_n = -1.0 * g.grad_nuc().reshape(self.mol.ndims) * self.FORCE_UNITS
            result['electronic_forces'] = f_e.defunits()
            result['nuclear_forces'] = f_n.defunits()
            result['forces'] = result['electronic_forces'] + result['nuclear_forces']

        # Return the result
        result['electronic_state'] = el_state
        self.last_el_state = el_state
        result['mulliken'] = DotDict({a: p for a, p in zip(self.mol.atoms, atom_pop)})
        result['mulliken'].type = 'atomic'
        return result

    def prep(self, force=False):
        # TODO: spin, isotopic mass, symmetry
        if self._prepped and not force: return
        self.pyscfmol = self._build_mol()
        self.kernel, self.reference = self._build_theories()
        self._prepped = True

    def _build_mol(self):
        """TODO: where does charge go? Model or molecule?"""
        pyscfmol = gto.Mole()
        pyscfmol.atom = [[atom.elem, atom.position.value_in(u.angstrom)]
                         for atom in self.mol.atoms]
        pyscfmol.basis = self._get_basis_name()
        pyscfmol.charge = self.get_formal_charge()
        if 'symmetry' in self.params and self.params.symmetry is not None:
            pyscfmol.symmetry = self.params.symmetry
        pyscfmol.stdout = self.logs
        with redirect_stderr(StringIO()) as builderr:
            pyscfmol.build()
        builderr.seek(0)
        for line in builderr:
            if line.strip() == 'Warn: Ipython shell catchs sys.args':
                continue
            else:
                self.logger.warning('PYSCF: ' + line)
        return pyscfmol

    def _converge(self, method, dm0=None):
        """
        Automatically try a bunch of fallback methods for convergence
        see also https://www.molpro.net/info/2015.1/doc/manual/node176.html#sec:difficulthf
        """
        # TODO: make this user configurable
        # TODO: generalize outside of pyscf

        energy = method.kernel(dm0=dm0)
        failed = []

        # stop here if it converged
        if method.converged:
            return method, failed

        # fallback 1: don't use previous density matrix OR change initial_guess
        failed.append(method)
        if dm0 is not None:
            method.init_guess = 'atom'
        else:
            method.init_guess = 'minao'
        self.logger.handled('SCF failed to converge. Retrying with initial guess %s' % method.init_guess)
        energy = method.kernel()
        if method.converged:
            return method, failed

        # fallback 2: level shift, slower convergence
        # this probably won't converge, but is intended to get us in the right basin for the next step
        # NEWFEATURE: should dynamically adjust level shift instead of hardcoded cycles
        self.logger.handled('SCF failed to converge. Performing %d iterations with level shift of -0.5 hartree'
                            % (method.max_cycle / 2))
        failed.append(method)
        method.init_guess = 'minao'
        method.level_shift = -0.5
        method.max_cycle /= 2
        energy = method.kernel()
        if method.converged:
            return method, failed

        # fallback 2 cont.: remove level shift and try to converge
        self.logger.handled('Removing level shift and continuing')
        level_shift_dm = method.make_rdm1()
        method.level_shift = 0.0
        method.max_cycle *= 2
        energy = method.kernel(dm0=level_shift_dm)
        if method.converged:
            return method, failed

        raise orbitals.ConvergenceError(method)

    def _build_theories(self):
        if self.params.theory in self.NEEDS_REFERENCE:
            scf_ref = if_not_none(self.params.scf_reference, 'rhf')
            reference = self.THEORIES[scf_ref](self.pyscfmol)
            if 'scf_cycles' in self.params:
                reference.max_cycle = self.params.scf_cycles
            reference.callback = StatusLogger('%s/%s reference procedure:' % (scf_ref, self.params.basis),
                                              ['cycle', 'e_tot'], self.logger)
        else:
            reference = None

        theory = self.THEORIES[self.params.theory](self.pyscfmol)
        theory.callback = StatusLogger('%s/%s procedure:' % (self.params.theory, self.params.basis),
                                          ['cycle', 'e_tot'], self.logger)

        if 'scf_cycles' in self.params:
            theory.max_cycle = self.params.scf_cycles
        return theory, reference

    def _get_ao_basis_functions(self):
        bfs = []
        for label in self.pyscfmol.cart_labels():
            atomidx, elem, shell, cart = label
            n = int(shell[0])
            l = orbitals.ANGMOM[shell[1:]]
            atom = self.mol.atoms[atomidx]
            assert atom.elem == elem
            bfs.append(orbitals.AtomicBasisFunction(atom, n=n, l=l, cart=cart))
        return bfs

    def _get_basis_name(self):
        """
        Translate basis_orbitals set name into a spec that pyscf recognizes
        :return:
        """
        # TODO: actually implement this
        return self.params.basis

    @staticmethod
    def _get_ao_matrices(mf):
        h1e = mf.get_hcore() * u.hartree
        sao = mf.get_ovlp()
        return DotDict(h1e=h1e, sao=sao)

    def _get_scf_matrices(self, mf, ao_mats):
        # TODO: density matrix units??? (dimensionless in this wfn, I think ...)
        dm = mf.make_rdm1()
        veff = mf.get_veff(dm=dm) * u.hartree
        fock = ao_mats.h1e + veff
        scf_matrices = dict(density_matrix=dm,
                            h2e=veff,
                            fock_ao=fock)
        scf_matrices.update(ao_mats)
        return scf_matrices


class StatusLogger(object):
    LEN = 15

    def __init__(self, description, columns, logger):
        self.logger = logger
        self.description = description
        self.columns = columns
        self._init = False
        self._row_format = ("{:<%d}" % self.LEN) + ("{:>%d}" % self.LEN) * (len(columns) - 1)

    def __call__(self, info):
        if not self._init:
            self.logger.status('Starting energy model calculation: %s' % self.description)
            self.logger.status(self._row_format.format(*self.columns))
            self.logger.status(self._row_format.format(*['-' * (self.LEN - 2) for i in self.columns]))
            self._init = True
        self.logger.status(self._row_format.format(*[info.get(c, 'n/a') for c in self.columns]))


class PySCFAOBasis(orbitals.AOBasis):

    def __getstate__(self):
        newdict = self.__dict__.copy()
        pmol = newdict.pop('pyscfmol', None)
        if pmol is not None:
            newdict['_pyscf_spec'] = dict(atom=pmol.atom, basis=pmol.basis)  # this is only for recreating the basis
        return newdict

    def __setstate__(self, state):
        spec = state.pop('_pyscf_spec')
        if spec is not None:
            state['pyscfmol'] = pyscf.gto.M(**spec)
        self.__dict__.update(state)

    def calculate_orb_grid(self, mo_in_ao, **kwargs):
        """
        Create and Populate a volumetric grid with the wfn values at each point
        :param mo_in_ao: 1-D array of molecular orbital coefficients in ao basis_orbitals
        :param kwargs: arguments for VolumetricGrid
        """
        # TODO: spatial wavefunction units??? length^(1/2), I think?
        spec = self.pyscfmol
        atom_positions = np.array([spec.atom_coord(i) for i in xrange(spec.natm)]) * u.angstrom
        grid = orbitals.VolumetricGrid(atom_positions, **kwargs)
        points = grid.xyzlist().reshape(3, grid.npoints ** 3).T.value_in(u.bohr)
        aovals = numint.eval_ao(spec, np.ascontiguousarray(points))
        movals = aovals.dot(mo_in_ao)
        grid.fxyz = movals.reshape(grid.npoints, grid.npoints, grid.npoints)
        return grid

    def get_eris_in_basis(self, orbs):
        """
        return electron repulsion integrals transformed in (in form eri[i,j,k,l] = (ij|kl))
        :param orbs:
        :return:
        """
        eri = ao2mo.full(self.pyscfmol, orbs.T, compact=True) * u.hartree
        eri.defunits_inplace()
        return ERI4FoldTensor(eri, orbs)


class ERI4FoldTensor(object):
    def __init__(self, mat, basis_orbitals):
        self.mat = mat
        self.basis_orbitals = basis_orbitals
        nbasis = len(self.basis_orbitals)
        mapping = np.zeros((nbasis, nbasis), dtype='int')
        ij = 0
        for i in xrange(nbasis):
            for j in xrange(i + 1):
                mapping[i, j] = mapping[j, i] = ij
                ij += 1
        self.mapping = mapping

    def __getitem__(self, item):
        i, j, k, l = item
        return self.mat[self.mapping[i, j], self.mapping[k, l]]
