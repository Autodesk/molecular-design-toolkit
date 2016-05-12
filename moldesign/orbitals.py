# NEW FEATURE: less nesting, more referencing
import collections
import ipywidgets as ipy
import numpy as np

from moldesign import units as u, ui
from moldesign.utils import if_not_none, DotDict, Alias
from moldesign.viewer import GeometryViewer


SHELLS = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h'}
ANGMOM = {v: k for k, v in SHELLS.iteritems()}
SPHERICALNAMES = {(0, 0): 's', (1, -1): 'p(x)', (1, 0): 'p(z)', (1, 1): 'p(y)',
                  (2, 0): 'd(z2)', (2, -2): 'd(xy)', (2, -1): 'd(yz)', (2, 1): 'd(xz)', (2, 2): 'd(x2-y2)'}


class ConvergenceError(Exception): pass


class AtomicBasisFunction(object):
    def __init__(self, atom, n=None, l=None, m=None, cart=None):
        self.atom = atom
        self.n = n
        self.l = l
        self.m = m
        if cart is not None:
            assert self.m is None, 'Both cartesian and spherical components passed!'
            assert len(cart) == self.l, 'Angular momentum does not match specified component %s' % cart
            for e in cart: assert e in 'xyz'
            self.cart = ''.join(sorted(cart))

    @property
    def orbtype(self):
        if self.l == 0: t = 's'
        elif self.cart is not None: t = SHELLS[self.l] + self.cart
        else: t = SPHERICALNAMES[self.l, self.m]
        return t

    @property
    def aotype(self):
        t = self.orbtype
        if self.n: return '%s%s' % (self.n, t)
        else: return t

    def __str__(self):
        return 'AO ' + self.name

    @property
    def name(self):
        try: return '%s on atom %s' % (self.aotype, self.atom.name)
        except: return 'Basis Fn'

    def __repr__(self):
        return '<%s %s>' % (self.__class__.__name__, self.name)


class AOBasis(object):
    """
    Stores an AO wfn description.
    RIGHT NOW, this has limited functionality, and should be subclassed to implement anything.
    """
    __len__ = Alias('basis_fns.__len__')
    __getindex__ = Alias('basis_fns.__getindex__')
    __setindex__ = Alias('basis_fns.__setindex__')
    __getitem__ = Alias('basis_fns.__getitem__')
    __setitem__ = Alias('basis_fns.__setitem__')

    def __init__(self, mol, basis_fns=None, name=None, h1e=None, overlaps=None,
                 **kwargs):
        self.mol = mol
        self.name = name
        self.h1e = h1e
        self.overlaps = overlaps
        self.basis_fns = basis_fns
        self.nbasis = len(self.basis_fns)
        for kw, val in kwargs.iteritems():
            setattr(self, kw, val)

        self.on_atom = {}
        if self.basis_fns is not None:
            for fn in self.basis_fns:
                self.on_atom.setdefault(fn.atom, []).append(fn)


    def calculate_orb_grid(self, coeffs, padding=2.5*u.angstrom, npoints=30):
        """Abstract method - subclasses need to implement this"""
        raise NotImplementedError('WFN grids are not implemented for this wfn set')

    def __repr__(self):
        return '<%s (%s) of %s>' % (self.__class__.__name__, self.name, self.mol)

    def as_mos(self, wfn=None):
        if wfn is not None: assert wfn.aobasis == self
        coeffs = np.identity(self.nbasis)
        orbitals = []
        for bsfn, c in zip(self.basis_fns, coeffs):
            orb = Orbital(c,
                          aobasis=self,
                          name=bsfn.name)
            orbitals.append(orb)
        atomic_orbs = MolecularOrbitals(orbitals, wfn=wfn, orbtype='AO basis')
        return atomic_orbs


class Orbital(object):
    """
    Stores a single orbital and its meta-data
    Generally wants to be part of a set of MolecularOrbitals
    """
    def __init__(self, coeffs,
                 wfn=None, aobasis=None, energy=None,
                 occupation=None, name='unnamed',
                 **kwargs):
        self.coeffs = np.array(coeffs)
        self.name = name
        self.energy = energy
        self.parent = None
        self.aobasis = aobasis
        self.occupation = occupation
        self.wfn = wfn
        if wfn is not None and wfn.aobasis is not None:
            if self.aobasis is not None:
                assert wfn.aobasis == self.aobasis
            else:
                self.aobasis = wfn.aobasis
        if self.aobasis is not None:
            assert len(self.coeffs) == len(self.aobasis)
        for k, v in kwargs.iteritems():
            setattr(self, k, v)

    def overlap(self, other):
        return self.dot(self.aobasis.overlaps.dot(other))

    def __repr__(self):
        return '<%s %s>' % (self.__class__.__name__, self.name)

    # This allows it to behave like a normal numpy array
    # Ideally we'd just subclass numpy.ndarray ... except that doesn't work
    __getitem__ = Alias('coeffs.__getitem__')
    __len__ = Alias('coeffs.__len__')
    shape = Alias('coeffs.shape')

    def __dir__(self):
        attributes = set(self.__dict__.keys() + dir(self.__class__) + dir(self.coeffs))
        return list(attributes)

    def __getattr__(self, item):
        if hasattr(self, 'coeffs'):
            return getattr(self.coeffs, item)
        else:
            raise AttributeError()


class MolecularOrbitals(object):
    """
    Stores a wfn of molecular orbitals in an AO wfn
    Orbitals are accessed as orbs[orbital index, ao index]
    """
    def __init__(self, orbitals, wfn=None,
                 canonical=False,
                 orbtype=None,
                 **kwargs):
        """
        :param coeffs: orbital coefficients (indexed as coeffs[ orbital_index, AO_basis_index])
        :param wfn: AOBasis-like object
        :param energies: energies of the orbitals
        :param orbtype: name of these orbitals (e.g. 'canonical' - the default, 'natural', 'nbo', etc.)
        """
        if canonical:
            assert orbtype is None or orbtype == 'canonical'
            orbtype = 'canonical'
        if orbtype == 'canonical': canonical = True
        self.orbitals = orbitals
        self.coeffs = np.array([orb.coeffs for orb in self.orbitals])
        self.wfn = wfn
        if self.wfn is None:
            self.wfn = self.orbitals[0].wfn
        self.orbtype = orbtype
        self.num_orbs = len(self.orbitals)
        if canonical: self._set_cmo_names()

        energies = self.fock.diagonal()
        for idx, orb in enumerate(self.orbitals):
            if orb.energy is None:
                orb.energy = energies[idx]
            else:
                # This is a loose check to make sure we're in the right ballpark
                if not np.allclose(orb.energy.defunits_value(),
                                   energies[idx].defunits_value(),
                                   atol=4):
                    raise ConvergenceError(
                            "Orbitals are not self-consistent; the SCF procedure likely did not converge")
            if orb.name is None:
                orb.name = '%s orb %d' % (self.orbtype, idx)

        for k, v in kwargs.iteritems():
            setattr(self, k, v)

    def align_phases(self, other, threshold = 0.5, assert_same_type=True):
        """
        Align the phases of these orbitals with those of the other
        :type other: MolecularOrbitals
        :param assert_same_type: require that both orbitals sets have the same orbtype
        """
        # TODO: we're using the wrong overlap matrix here, since the atoms have moved. Need to get overlaps
        # between the two steps for this to work consistently
        if assert_same_type:
            assert self.orbtype == other.orbtype, "Orbital type mismatch: %s vs. %s" % (self.orbtype, other.orbtype)
        for thisorb, otherorb in zip(self, other):
            olap = thisorb.overlap(otherorb)
            if thisorb.overlap(otherorb) < -1.0 * threshold:
                thisorb.coeffs *= -1.0
                # TODO: print a warning if overlap is small?

    def overlap(self, other):
        return self.dot(self.aobasis.overlaps.dot(other.T))

    def calculate_orb_grid(self, orbidx, **kwargs):
        """
        :param orbidx: orbital index
        :param kwargs: arguments to AOBasis
        :return type: VolumetricGrid
        """
        return self.wfn.aobasis.calculate_orb_grid(self[orbidx], **kwargs)

    def __iter__(self):
        return iter(self.orbitals)

    def __str__(self):
        return '%s orbitals' % self.orbtype

    def __repr__(self):
        return '<%d %s %s in %s>' % (self.num_orbs, self.orbtype,
                                     self.__class__.__name__, str(self.wfn))

    @property
    def fock(self):
        return self.from_ao(self.wfn.fock_ao)

    @property
    def overlaps(self):
        return self.from_ao(self.wfn.aobasis.overlaps)

    @property
    def h1e(self):
        return self.from_ao(self.wfn.aobasis.h1e)

    @property
    def h2e(self):
        return self.fock - self.h1e

    def from_ao(self, ao_operator):
        """
        Transform an operator from the AO basis
        Dot doesn't place nice with units, so we need to pass them explicitly
        """
        ao_units = u.get_units(ao_operator)
        return self.coeffs.dot(ao_operator.dot(self.coeffs.T)) * ao_units

    def to_ao(self, mo_operator):
        """
        Rotate an operator from mo basis to (generally non-orthogonal) ao basis
        Dot doesn't place nice with units, so we need to pass them explicitly
        """
        units = u.get_units(mo_operator)
        s = self.wfn.aobasis.overlaps
        o_ao = s.dot(self.coeffs.T).dot(mo_operator).dot(self.coeffs).dot(s)
        return o_ao * units

    def _set_cmo_names(self):
        for i, orb in enumerate(self.orbitals):
            if orb.name != 'unnamed' and orb.name is not None:
                continue
            if i <= self.wfn.homo:
                if i < self.wfn.homo - 2:
                    orb.name = 'cmo %d' % i
                elif i == self.wfn.homo:
                    orb.name = 'HOMO'
                else:
                    orb.name = 'HOMO-%d' % (self.wfn.homo - i)
            else:
                if i == self.wfn.lumo:
                    orb.name = 'LUMO'
                elif i <= self.wfn.lumo + 2:
                    orb.name = 'LUMO+%d' % (i - self.wfn.lumo)
                else:
                    orb.name = 'virt cmo %d' % i

    # This allows it to behave like a normal numpy array
    # Ideally we'd just subclass numpy.ndarray ... except that doesnt work
    __getitem__ = Alias('orbitals.__getitem__')
    __len__ = Alias('orbitals.__len__')
    shape = Alias('coeffs.shape')

    def __dir__(self):
        attributes = set(self.__dict__.keys() + dir(self.__class__) +
                         dir(self.coeffs) + dir(self.orbitals[0]))
        return list(attributes)

    def __getattr__(self, item):
        if 'coeffs' in self.__dict__:
            try: return getattr(self.coeffs, item)
            except AttributeError: pass

        if 'orbitals' in self.__dict__:
            try:
                orbprops = [getattr(orb, item) for orb in self.orbitals]
            except AttributeError:
                raise AttributeError('%s has no attribute %s' % (self, item))

            try: return u.to_units_array(orbprops)
            except TypeError: return orbprops

        else:
            raise AttributeError()


class ElectronicState(object):
    """
    Stores the results of a quantum chemistry calculation.
    This is necessarily pretty flexible, but generally stores an AO wfn and one or more sets of orbitals.
    Can also store CI vectors, etc.
    """

    def __init__(self, mol, num_electrons,
                 theory=None,
                 aobasis=None, fock_ao=None,
                 civectors=None, **kwargs):
        """
        :param mol:
        :param aobasis:
        :param canonical_orbitals:
        :param kwargs:
        :return:
        """
        self.mol = mol
        self.theory = theory
        self.civectors = civectors
        self.aobasis = aobasis
        self.orbitals = DotDict()
        self.fock_ao = fock_ao
        self.num_electrons = num_electrons
        self.homo = self.num_electrons/2 - 1
        self.lumo = self.homo + 1
        self._has_canonical = False

        if self.aobasis is not None:
            self.orbitals['atomic'] = self.aobasis.as_mos(wfn=self)

        for arg, val in kwargs.iteritems():
            setattr(self, arg, val)

    def __repr__(self):
        return '<ElectronicState (%s) of %s>' % (self.description, str(self.mol))

    def __str__(self):
        return '%s wfn' % self.description

    @property
    def description(self):
        return '%s/%s' % (self.theory, self.aobasis.name)

    def set_canonical_mos(self, orbs):
        if orbs.wfn is None: orbs.wfn = self
        if self.fock_ao is None and orbs.energy is not None:
            fock_cmo = orbs.energy * np.identity(len(self.aobasis))
            self.fock_ao = orbs.to_ao(fock_cmo)
        self._has_canonical = True

    def align_orbital_phases(self, other, assert_same=True):
        """Align this wavefunction's orbitals to have the same phase as those in `other`.
        :type other: ElectronicState
        :param assert_same: raise an exception if the two wavefunctions do not have the same kinds of orbitals
        """
        for orbtype in self.orbitals:
            if orbtype not in other.orbitals:
                if assert_same: assert False, '%s has orbital type %s, but %s does not.' % (self, orbtype, other)
                else: continue
            self.orbitals[orbtype].align_phases(other.orbitals[orbtype])

    def run_nbo(self, **kwargs):
        from moldesign.interfaces import nbo_interface
        nbo_interface.run_nbo(self.mol, **kwargs)

    def add_orbitals(self, orbs, orbtype='canonical', **kwargs):
        mo_object = MolecularOrbitals(orbs,
                                      wfn=self,
                                      orbtype=orbtype)
        self.orbitals[orbtype] = mo_object
        if orbtype == 'canonical' and not self._has_canonical:
            self.set_canonical_mos(mo_object)
        return mo_object

    @property
    def molecular_orbitals(self):
        """A synonym for self.orbitals['canonical'], since this is usually what's wanted"""
        return self.orbitals['canonical']

    @molecular_orbitals.setter
    def molecular_orbitals(self, val):
        """A synonym for self.orbitals['canonical'], since this is usually what's wanted"""
        self.orbitals['canonical'] = val


class VolumetricGrid(object):
    """
    Helper object for preparing gaussian CUBE files
    """
    UNITS = u.angstrom
    def __init__(self, positions, padding=2.5*u.angstrom, npoints=25):
        mins = positions.min(axis=0) - padding
        maxes = positions.max(axis=0) + padding
        self.npoints = npoints
        self.xr = (mins[0], maxes[0])
        self.yr = (mins[1], maxes[1])
        self.zr = (mins[2], maxes[2])
        self._origin = mins.value_in(self.UNITS)
        self.dx = (self.xr[1] - self.xr[0]).value_in(self.UNITS) / (float(npoints) - 1)
        self.dy = (self.yr[1] - self.yr[0]).value_in(self.UNITS) / (float(npoints) - 1)
        self.dz = (self.zr[1] - self.zr[0]).value_in(self.UNITS) / (float(npoints) - 1)
        self.fxyz = None

    def xyzlist(self):
        stride = self.npoints * 1j
        grids = np.mgrid[self.xr[0]:self.xr[1]:stride,
                self.yr[0]:self.yr[1]:stride,
                self.zr[0]:self.zr[1]:stride]
        return grids * self.UNITS

    def origin(self):
        return tuple(self._origin)


class OrbitalViz(ui.SelectionGroup):
    def __init__(self, mol, **kwargs):
        """
        :param mol: a molecule with A) orbitals, and B) an energy model with calculate_orbital_grid
        :param kwargs: kwargs for the viewer
        :return:
        """
        self.viewer = GeometryViewer(mol=mol, **kwargs)
        self.viewer.wfn = mol.electronic_state
        self.uipane = OrbitalUIPane(self, height=int(self.viewer.height)-50)
        hb = ipy.HBox([self.viewer, self.uipane])
        super(OrbitalViz, self).__init__([hb])


class OrbitalUIPane(ui.Selector, ipy.Box):
    # TODO: deal with orbitals not present in all frames of a trajectory
    # TODO: deal with orbital properties changing over a trajectory
    def __init__(self, viz, **kwargs):
        self.viz = viz
        kwargs.setdefault('width', 325)

        self.type_dropdown = ipy.Dropdown(options=self.viz.viewer.wfn.orbitals.keys())
        initial_orb = 'canonical'
        if initial_orb not in self.type_dropdown.options:
            initial_orb = self.type_dropdown.options.iterkeys().next()
        self.type_dropdown.value = initial_orb
        self.type_dropdown.observe(self.new_orb_type, 'value')

        self.orblist = ipy.Select(width=kwargs['width'], height=int(kwargs['height']) - 75)

        self.isoval_selector = ui.create_value_selector(ipy.FloatSlider,
                                                        value_selects='orbital_isovalue',
                                                        min=0.0, max=0.05,
                                                        value=0.01, step=0.0001,
                                                        width=kwargs['width'],
                                                        description='Isovalue')

        self.orb_resolution = ipy.Text(description='Orbital resolution', width=75)
        self.orb_resolution.value = '40'  # this is a string to enable the 'on_submit' method
        self.change_resolution()
        self.orb_resolution.on_submit(self.change_resolution)

        children = [self.type_dropdown, self.orblist, self.isoval_selector, self.orb_resolution]
        super(OrbitalUIPane, self).__init__(children, **kwargs)
        self.new_orb_type()
        self.orblist.observe(self.new_orbital_selection, 'value')


    def new_orbital_selection(self, *args):
        self.fire_selection_event({'orbname': (self.type_dropdown.value, self.orblist.value)})

    def handle_selection_event(self, *args):
        # TODO: update the selected orbitals if something actually else triggers this
        pass

    def new_orb_type(self, *args):
        """Create list of available orbitals when user selects a new type
        """
        wfn = self.viz.viewer.wfn
        newtype = self.type_dropdown.value
        neworbs = wfn.orbitals[newtype]
        orblist = collections.OrderedDict()

        for i, orb in enumerate(neworbs):
            if hasattr(orb, 'unicode_name'):
                orbname = orb.unicode_name
            else:
                orbname = orb.name

            meta = ''
            if orb.energy is not None:
                meta = '{:.02fP}'.format(orb.energy.defunits())
            if orb.occupation is not None:
                if meta: meta += ', '
                meta += 'occ %.2f' % orb.occupation
            if meta:
                desc = '%d. %s   (%s)' % (i, orbname, meta)
            else:
                desc = '%d. %s' % (i, orbname)
            orblist[desc] = i
        self.orblist.options = orblist

    def change_resolution(self, *args):
        viewer = self.viz.viewer
        viewer._orbital_kwargs['npts'] = int(self.orb_resolution.value)
        if viewer.current_orbital is not None:
            viewer.draw_orbital(viewer.current_orbital, render=True, **viewer._orbital_kwargs)
