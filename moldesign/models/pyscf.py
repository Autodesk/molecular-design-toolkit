# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from __future__ import absolute_import  # prevent clashes between this and the "pyscf" package
from cStringIO import StringIO

import numpy as np

import moldesign as mdt
from moldesign import units as u
from moldesign import compute, orbitals
from moldesign.interfaces.pyscf_interface import force_remote, mol_to_pyscf, \
    StatusLogger, SPHERICAL_NAMES
from .base import QMBase
from moldesign import uibase
from moldesign.utils import DotDict


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


class LazyClassMap(object):
    """ For lazily importing classes from modules (when there's a lot of import overhead)

    Class names should be stored as their *absolute import strings* so that they can be imported
    only when needed

    Example:
        >>> myclasses = LazyClassMap({'od': 'collections.OrderedDict'})
        >>> myclasss['od']()
        OrderedDict()
    """
    def __init__(self, mapping):
        self.mapping = mapping

    def __getitem__(self, key):
        import importlib
        fields = self.mapping[key].split('.')
        cls = fields[-1]
        modname = '.'.join(fields[:-1])
        mod = importlib.import_module(modname)
        return getattr(mod, cls)

    def __contains__(self, item):
        return item in self.mapping

    def __iter__(self):
        return iter(self.mapping)

# PySCF metadata constants
THEORIES = LazyClassMap({'hf': 'pyscf.scf.RHF', 'rhf': 'pyscf.scf.RHF',
                         'uhf': 'pyscf.scf.UHF',
                         'mcscf': 'pyscf.mcscf.CASSCF', 'casscf': 'pyscf.mcscf.CASSCF',
                         'casci': 'pyscf.mcscf.CASCI',
                         'mp2': 'pyscf.mp.MP2',
                         'dft': 'pyscf.dft.RKS', 'rks': 'pyscf.dft.RKS', 'ks': 'pyscf.dft.RKS'})

NEEDS_REFERENCE = set('mcscf casscf casci mp2'.split())
NEEDS_FUNCTIONAL = set('dft rks ks uks'.split())
IS_SCF = set('rhf uhf hf dft rks ks'.split())
FORCE_CALCULATORS = LazyClassMap({'rhf': 'pyscf.grad.RHF', 'hf': 'pyscf.grad.RHF',
                                  'rks': 'pyscf.grad.RKS', 'ks': 'pyscf.grad.RKS'})


@exports
class PySCFPotential(QMBase):
    DEFAULT_PROPERTIES = ['potential_energy',
                          'wfn',
                          'mulliken']
    ALL_PROPERTIES = DEFAULT_PROPERTIES + ['eri_tensor',
                                           'forces',
                                           'nuclear_forces',
                                           'electronic_forces']
    PARAM_SUPPORT = {'theory': ['rhf', 'rks', 'mp2'],
                     'functional': ['b3lyp', 'blyp', 'pbe0', 'x3lyp', 'MPW3LYP5']}

    FORCE_UNITS = u.hartree / u.bohr

    @mdt.utils.kwargs_from(QMBase)
    def __init__(self, **kwargs):
        super(PySCFPotential, self).__init__(**kwargs)
        self.pyscfmol = None
        self.reference = None
        self.kernel = None
        self.logs = StringIO()
        self.logger = uibase.Logger('PySCF interface')

    @compute.runsremotely(enable=force_remote, is_imethod=True)
    def calculate(self, requests=None):
        self.logger = uibase.Logger('PySCF calc')
        do_forces = 'forces' in requests
        if do_forces and self.params.theory not in FORCE_CALCULATORS:
            raise ValueError('Forces are only available for the following theories:'
                             ','.join(FORCE_CALCULATORS))
        if do_forces:
            force_calculator = FORCE_CALCULATORS[self.params.theory]

        self.prep(force=True)  # rebuild every time

        # Set up initial guess
        if self.params.wfn_guess == 'stored':
            dm0 = self.params.initial_guess.density_matrix_ao
        else:
            dm0 = None

        # Compute reference WFN (if needed)
        refobj = self.pyscfmol
        if self.params.theory in NEEDS_REFERENCE:
            reference = self._build_theory(self.params.get('reference', 'rhf'),
                                           refobj)
            kernel, failures = self._converge(reference, dm0=dm0)
            refobj = self.reference = kernel
        else:
            self.reference = None

        # Compute WFN
        theory = self._build_theory(self.params['theory'],
                                    refobj)
        if self.params['theory'] not in IS_SCF:
            theory.kernel()
            self.kernel = theory
        else:
            self.kernel, failures = self._converge(theory, dm0=dm0)

        # Compute forces (if requested)
        if do_forces:
            grad = force_calculator(self.kernel)
        else:
            grad = None

        props = self._get_properties(self.reference, self.kernel, grad)

        if self.params.store_orb_guesses:
            self.params.wfn_guess = 'stored'
            self.params.initial_guess = props['wfn']

        return props

    def _get_properties(self, ref, kernel, grad):
        """ Analyze calculation results and return molecular properties

        Args:
            ref (pyscf.Kernel): Reference kernel (can be None)
            kernel (pyscf.Kernel): Theory kernel
            grad (pyscf.Gradient): Gradient calculation

        Returns:
            dict: Molecular property names and values
        """
        result = {}

        if self.reference is not None:
            result['reference_energy'] = (ref.e_tot*u.hartree).defunits()
            # TODO: check sign on correlation energy. Is this true for anything besides MP2?
            if hasattr(kernel, 'e_corr'):
                result['correlation_energy'] = (kernel.e_corr *u.hartree).defunits()
                result['potential_energy'] = result['correlation_energy'] +\
                                             result['reference_energy']
            else:
                result['potential_energy'] = (kernel.e_tot*u.hartree).defunits()

            orb_calc = ref
        else:
            result['potential_energy'] = (kernel.e_tot*u.hartree).defunits()
            orb_calc = kernel

        if grad is not None:
            f_e = -1.0 * grad.grad_elec() * self.FORCE_UNITS
            f_n = -1.0 * grad.grad_nuc() * self.FORCE_UNITS
            result['electronic_forces'] = f_e.defunits()
            result['nuclear_forces'] = f_n.defunits()
            result['forces'] = result['electronic_forces'] + result['nuclear_forces']

        if self.params.theory in ('casscf', 'mcscf'):
            from pyscf.mcscf import CASCI
            casobj = CASCI(ref,
                           self.params.active_orbitals,
                           self.params.active_electrons)
        elif self.params.theory == 'casci':
            casobj = kernel

        if self.params.theory in ('casscf', 'mcscf', 'casci'):
            orb_calc = kernel  # get the orbs directly from the post-HF theory
            casobj.fcisolver.nroots = self.params.get('num_states',
                                                      self.params.state_average)
            casresult = casobj.kernel()
            result['state_energies'] = (casresult[0] * u.hartree).defunits()
            result['ci_vectors'] = map(self._parse_fci_vector, casresult[2])

            # potential_energy is the energy of the molecule's assigned state
            result['state_averaged_energy'] = result['potential_energy']
            result['potential_energy'] = result['state_energies'][self.mol.electronic_state_index]
            # TODO: add 'reference wavefunction' to result
            ao_obj = ref
            dips, tdips = _get_multiconf_dipoles(self.pyscfmol, casobj, len(casobj.ci))
            result['state_dipole_moments'] = dips
            result['transition_dipole_moments'] = tdips
            result['dipole_moment'] = dips[0]

            # TODO: this is general, put it somewhere else
            oscs = {}
            nstates = len(result['state_energies'])
            for i in xrange(nstates):
                for j in xrange(i+1, nstates):
                    excitation_energy = result['state_energies'][j]-result['state_energies'][i]
                    tdip = result['transition_dipole_moments'][i, j].norm()
                    oscs[i, j] = (2.0*tdip ** 2*u.m_e*excitation_energy/
                                  (3.0*u.q_e ** 2*u.hbar ** 2)).to(u.ureg.dimensionless).magnitude
                    oscs[j, i] = -oscs[i, j]
            result['oscillator_strengths'] = oscs
        else:
            ao_obj = orb_calc

        ao_matrices = self._get_ao_matrices(ao_obj)
        scf_matrices = self._get_scf_matrices(orb_calc, ao_matrices)
        if hasattr(orb_calc, 'mulliken_pop'):
            ao_pop, atom_pop = orb_calc.mulliken_pop(verbose=-1)
            result['mulliken'] = DotDict({a: p for a, p in zip(self.mol.atoms, atom_pop)})
            result['mulliken'].type = 'atomic'

        if hasattr(orb_calc, 'dip_moment'):
            result['dipole_moment'] = orb_calc.dip_moment() * u.debye

        # Build the electronic state object
        basis = orbitals.basis.BasisSet(self.mol,
                                        orbitals=self._get_ao_basis_functions(),
                                        h1e=ao_matrices.h1e.defunits(),
                                        overlaps=scf_matrices.pop('sao'),
                                        name=self.params.basis)
        el_state = orbitals.wfn.ElectronicWfn(self.mol,
                                              self.pyscfmol.nelectron,
                                              aobasis=basis,
                                              fock_ao=scf_matrices['fock_ao'],
                                              density_matrix_ao=scf_matrices['density_matrix_ao'],
                                              description=self.theoryname)

        # Build and store the canonical orbitals
        cmos = []
        for iorb, (coeffs, energy) in enumerate(zip(orb_calc.mo_coeff.T,
                                                    orb_calc.mo_energy * u.hartree)):
            cmos.append(orbitals.Orbital(coeffs, wfn=el_state))
        if hasattr(orb_calc, 'get_occ'):
            for orb, occ in zip(cmos, orb_calc.get_occ()):
                orb.occupation = occ
        el_state.add_orbitals(cmos, orbtype='canonical')

        # Return the result
        result['wfn'] = el_state
        return result

    def prep(self, force=False):
        # TODO: spin, isotopic mass, symmetry
        for p in 'basis theory'.split():
            if self.params.get(p, None) is None:
                raise ValueError('Parameter "%s" is required' % p)

        if self._prepped and not force: return
        self.pyscfmol = self._build_mol()
        self._prepped = True

    def _build_mol(self):
        """TODO: where does charge go? Model or molecule?"""
        pyscfmol = mol_to_pyscf(self.mol, self.params.basis,
                                symmetry=self.params.get('symmetry', None),
                                charge=self.get_formal_charge())
        pyscfmol.stdout = self.logs
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

        raise mdt.QMConvergenceError(method)

    def _build_theory(self, name, refobj):
        if name in ('mscscf', 'casci', 'casscf'):
            theory = THEORIES[name](refobj,
                                    self.params.active_orbitals,
                                    self.params.active_electrons)

            if name != 'casci' and self.params.state_average > 1:
                theory = theory.state_average_([1.0/self.params.state_average
                                                for i in xrange(self.params.state_average)])
        else:
            theory = THEORIES[name](refobj)

        theory.callback = StatusLogger('%s procedure:' % self.theoryname,
                                       ['cycle', 'e_tot'],
                                       self.logger)

        if 'scf_cycles' in self.params:
            theory.max_cycle = self.params.scf_cycles

        if 'functional' in self.params:
            self._assign_functional(theory, name,
                                    self.params.get('functional', None))

        return theory

    _OCCMAP = {('0', '0'): '0',
               ('1', '0'): 'a',
               ('0', '1'): 'b',
               ('1', '1'): '2'}

    @property
    def theoryname(self):
        p = self.params
        if p.theory == 'rks':
            th = 'RKS(%s)' % p.functional.upper()
        elif p.theory in ('casscf', 'casci'):
            th = '%s(%d,%d)' % (p.theory.upper(), p.active_orbitals, p.active_electrons)
            if p.theory == 'casscf' and p.state_average > 1:
                th += ' SA-%d' % p.state_average
        else:
            th = p.theory.upper()

        return '%s/%s' % (th, p.basis)

    def _parse_fci_vector(self, ci_vecmat):
        """ Translate the PySCF FCI matrix into a dictionary of configurations and weights

        Args:
            ci_vecmat (np.ndarray): ci vector from a PySCF FCI calculation

        Returns:
            Mapping[str, float]: dictionary of configuration weights (normalized) organized by
                configuration label. Configurations labeled by their active space orbital
                occupations: 0 (unoccupied), a (alpha electron only), b (beta electron only), or '2'
                (doubly occupied)

        Example:
            >>> import numpy as np
            >>> model = PySCFPotential(active_orbitals=2, active_electrons=2)
            >>> model._parse_fci_vector(np.array([[1.0, 2.0],[3.0, 4.0]]))
            {'20': 1.0,
             'ba': 2.0,
             'ab': 3.0,
             '02': 4.0}
        """
        from pyscf.fci import cistring
        conf_bin = cistring.gen_strings4orblist(range(self.params.active_orbitals),
                                                self.params.active_electrons/2)
        civecs = {}
        for i, ca in enumerate(conf_bin):
            for j, cb in enumerate(conf_bin):
                astring = bin(ca)[2:].zfill(self.params.active_orbitals)
                bstring = bin(cb)[2:].zfill(self.params.active_orbitals)
                s = ''.join(reversed([self._OCCMAP[a, b] for a, b in zip(astring, bstring)]))
                civecs[s] = ci_vecmat[i, j]
        return civecs

    @staticmethod
    def _assign_functional(kernel, theory, fname):
        if theory in NEEDS_FUNCTIONAL:
            if fname is not None:
                kernel.xc = fname
            else:
                raise ValueError('No functional specified for reference theory "%s"' % theory)

    def _get_ao_basis_functions(self):
        """ Convert pyscf basis functions into a list of atomic basis functions

        Notes:
            PySCF stores *shells* instead of a flat list, so we need to do a little hacky
                guesswork to do this conversion. We include consistentcy checks with the annotated
                list of basis functions stored from ``mole.cart_labels()``
            As of PySCF v1.0, only cartesian orbitals appear to be supported, and that's all
                supported here right now

        Returns:
            List[moldesign.Gaussians.AtomicBasisFunction]
        """
        bfs = []
        pmol = self.pyscfmol

        orblabels = iter(pmol.spheric_labels())

        for ishell in xrange(pmol.nbas):  # loop over shells (n,l)
            atom = self.mol.atoms[pmol.bas_atom(ishell)]
            angular = pmol.bas_angular(ishell)
            num_momentum_states = angular*2 + 1
            exps = pmol.bas_exp(ishell)
            num_contractions = pmol.bas_nctr(ishell)
            coeffs = pmol.bas_ctr_coeff(ishell)

            for ictr in xrange(num_contractions):  # loop over contractions in shell
                for ibas in xrange(num_momentum_states):  # loop over angular states in shell
                    label = orblabels.next()
                    sphere_label = label[3]
                    l, m = SPHERICAL_NAMES[sphere_label]
                    assert l == angular
                    # TODO: This is not really the principal quantum number
                    n = int(''.join(x for x in label[2] if x.isdigit()))

                    primitives = [orbitals.SphericalGaussian(atom.position.copy(),
                                                       exp, n, l, m,
                                                       coeff=coeff[ictr])
                                  for exp, coeff in zip(exps, coeffs)]
                    bfs.append(orbitals.AtomicBasisFunction(atom, n=n, l=angular, m=m,
                                                      primitives=primitives))

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
        dm = mf.make_rdm1()
        veff = mf.get_veff(dm=dm) * u.hartree
        fock = ao_mats.h1e + veff
        scf_matrices = dict(density_matrix_ao=dm,
                            h2e=veff,
                            fock_ao=fock)
        scf_matrices.update(ao_mats)
        return scf_matrices


def _get_multiconf_dipoles(basis, mcstate, nstates):
    """ Compute dipoles and transition dipoles. Adapted from PySCF examples

    Note:
        Dipole moments are computed using the center of the nuclear charge as the origin. Dipole
        moments will need to be annotated or translated appropriately for charges systems.

    Args:
        basis ():
        mcstate ():
        nstates ():

    Returns:
        List[u.Vector[dipole]]: Dipole moments for each state
        Mapping[Tuple[int, int], u.Vector[dipole]]: mapping from pairs of state ids to transition
           dipole moments

    References:
        https://github.com/sunqm/pyscf/blob/e4d824853c49b7c19eb35cd6f9fe6ea675de932d/examples/1-advanced/030-transition_dipole.py
    """
    nuc_charges = [basis.atom_charge(i) for i in xrange(basis.natm)]
    nuc_coords = [basis.atom_coord(i) for i in xrange(basis.natm)]
    nuc_charge_center = np.einsum('z,zx->x', nuc_charges, nuc_coords)/sum(nuc_charges)
    basis.set_common_orig_(nuc_charge_center)
    nuc_dip = np.einsum('i,ix->x', nuc_charges, nuc_coords-nuc_charge_center) * u.a0 * u.q_e

    dip_ints = basis.intor('cint1e_r_sph', comp=3)
    orbcas = mcstate.mo_coeff[:, mcstate.ncore:mcstate.ncore+mcstate.ncas]

    dipoles, transitions = [], {}
    for istate in xrange(nstates):
        for fstate in xrange(istate, nstates):
            t_dm1 = mcstate.fcisolver.trans_rdm1(mcstate.ci[istate], mcstate.ci[fstate],
                                                 mcstate.ncas, mcstate.nelecas)
            t_dm1_ao = reduce(np.dot, (orbcas, t_dm1, orbcas.T))
            moment = np.einsum('xij,ji->x', dip_ints, t_dm1_ao) * u.a0 * u.q_e

            if istate == fstate:
                dipoles.append(moment)
            else:
                transitions[istate, fstate] = transitions[fstate, istate] = moment

    for idip, d in enumerate(dipoles):
        dipoles[idip] = nuc_dip - d

    return dipoles, transitions
