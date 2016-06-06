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

""" Class definitions for atomic and molecular orbitals.

Notes:
    In this documentation, we use the following conventions for labeling orbitals:
      - atomic orbitals using lower case greek labels and subscripts, e.g.,
          :math:`\left| \mu \right \rangle, F_{\nu \lambda}, etc.
      - molecular orbitals use lower case labels and subscripts, e.g.,
          :math:`\left| i \right \rangle, F_{kl}, etc.
      - adiabatic electronic states are indexed using capital letters, _N_, _L_, _M_, etc.

"""
import collections

import ipywidgets as ipy
import numpy as np

from moldesign import units as u, ui
from moldesign.utils import DotDict, Alias
from moldesign.viewer import GeometryViewer

SHELLS = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h'}
ANGMOM = {v: k for k, v in SHELLS.iteritems()}
SPHERICALNAMES = {(0, 0): 's', (1, -1): 'p(x)', (1, 0): 'p(z)', (1, 1): 'p(y)',
                  (2, 0): 'd(z^2)', (2, -2): 'd(xy)', (2, -1): 'd(yz)', (2, 1): 'd(xz)', (2, 2): 'd(x^2-y^2)'}


class ConvergenceError(Exception): pass


class PrimitiveBase(object):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError()

    def __call__(self, position):
        raise NotImplementedError()

    def overlap(self, other):
        raise NotImplementedError()

    @property
    def norm(self, other):
        raise NotImplementedError()


class Orbital(object):
    r"""
    Stores a single orbital and its meta-data
    Generally wants to be part of a set of MolecularOrbitals

    The orbital is defined as
    .. math::
        \left| i \right \rangle = \sum_\mu c_{i \mu} \left| \mu \right \rangle

    where the coefficients :math:`c_{i \mu}` are stored in ``self.coeffs`` and the basis orbitals
    :math:`\left| \mu \right \rangle` are stored at ``self.basis``
    """
    def __init__(self, coeffs, basis=None, wfn=None,
                 occupation=None, name='unnamed',
                 **kwargs):
        """ Initialization:

        Args:
            coeffs (numpy.array): orbital coefficients
            basis (BasisSet): basis for this orbital
            wfn (ElectronicState): total electronic wavefunction that this orbital is a part of
            occupation (float): occupation of this orbital
            name (str): optional label for this orbital:
        """
        self.coeffs = np.array(coeffs)
        self.name = name
        self.parent = None
        self.basis = basis
        self.occupation = occupation
        self.wfn = wfn

        # Assign the basis functions
        if wfn is not None and wfn.aobasis is not None:
            if self.basis is not None:
                assert wfn.aobasis is self.basis
            else:
                self.basis = wfn.aobasis
        if self.basis is not None:
            assert len(self.coeffs) == len(self.basis)

        # Assign arbitrary attributes
        for k, v in kwargs.iteritems():
            setattr(self, k, v)

    def overlap(self, other):
        """ Calculate overlap with another orbital

        Args:
            other (Orbital): calculate the overlap with this orbital

        Returns:
            float: orbital overlap
        """
        return self.dot(self.basis.overlaps.dot(other))


    def fock_element(self, other):
        """ Calculate fock matrix element with another orbital

        Args:
            other (Orbital): calculate the fock element with this orbital

        Returns:
            u.Scalar[energy]: fock matrix element
        """
        return self.coeffs.dot(self.wfn.fock_ao.dot(other.coeffs))

    @property
    def energy(self):
        """ u.Scalar[energy]: This orbital's energy

        Note:
            This is equivalent to self.fock(self)
        """
        return self.fock(self)

    def __call__(self, coords):
        """ Calculate the orbital's value at a position in space

        Args:
            coords (u.Vector[length]): Coordinates

        Returns:
            u.Scalar[length**(-3/2)]: value of the orbital
        """
        val = 0.0
        for c, fn in zip(self.coeffs, self.basis.orbitals):
            val += c * fn(coords)
        return val

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
    def __init__(self, orbitals, wfn=None, basis=None,
                 canonical=False,
                 orbtype=None,
                 **kwargs):
        """ Initialization:

        Args:
            orbitals (List[Orbital] OR numpy.array): EITHER a list of orbitals OR an
               array of coefficients (indexed as coeffs[orbital_index, basis_index])
            basis (BasisSet): orbital basis set
            wfn (ElectronicState): total electronic wavefunction that this orbital is a part of
            canonical (bool): designate these as the canonical orbitals
            orbtype (str): description of these orbitals
        """
        # Determine if these are the canonical orbitals
        if canonical:
            assert orbtype is None or orbtype == 'canonical'
            orbtype = 'canonical'
        if orbtype == 'canonical':
            canonical = True

        if not hasattr(orbitals[0], 'basis'):
            coeffs = orbitals
            orbitals = []
            for c in coeffs:
                orbitals.append(Orbital(c, basis=basis, wfn=wfn))

        self.orbitals = orbitals
        self.coeffs = np.array([orb.coeffs for orb in self.orbitals])
        self.wfn = wfn
        self.basis = basis
        self.orbtype = orbtype

        if self.wfn is None:
            self.wfn = self.orbitals[0].wfn
        if self.basis is None:
            self.basis = self.orbitals[0].basis

        for iorb, orbital in enumerate(self.orbitals):
            orbital.index = iorb
            assert orbital.basis == self.basis
            assert orbital.wfn == self.wfn
            orbital.coeffs = self.coeffs[iorb, :]

        if canonical:
            self._set_cmo_names()

        self._check_orb_energies()

        for k, v in kwargs.iteritems():
            setattr(self, k, v)

    def _check_orb_energies(self):
        energies = self.energies

        for idx, orb in enumerate(self.orbitals):
            # This is a loose check to make sure we're in the right ballpark
            if not np.allclose(orb.energy.defunits_value(),
                               energies[idx].defunits_value(),
                               atol=4):
                raise ConvergenceError(
                        "Orbitals are not self-consistent; the SCF procedure likely did not converge")
        if orb.name is None:
            orb.name = '%s orb %d'%(self.orbtype, idx)

    def align_phases(self, other, threshold=0.5, assert_same_type=True):
        """ Flip the signs of these orbitals to bring them into maximum coincidence
        with another set of orbitals

        Args:
            other (MolecularOrbitals): the "reference" set of orbitals to match phases with
            threshold (float): only flip orbital if the overlap is less than -1*threshold
            assert_same_type (bool): require that ``self.orbtype == other.orbtype``

        Note:
            This function assumes that the overlap matrix is the same for both sets of
            orbitals - this is a reasonable assumption if the two sets of orbitals were
            calculated at very similar molecular geometries.
        """
        if assert_same_type:
            assert self.orbtype == other.orbtype, "Orbital type mismatch: %s vs. %s" % (
                self.orbtype, other.orbtype)
        for thisorb, otherorb in zip(self, other):
            olap = thisorb.overlap(otherorb)
            if thisorb.overlap(otherorb) < -1.0 * threshold:
                thisorb.coeffs *= -1.0
                # TODO: print a warning if overlap is small?

    def overlap(self, other):
        """ Calculate overlaps between this and another set of orbitals

        Args:
            other (MolecularOrbitals):

        Returns:
            numpy.array: overlaps between the two sets of orbitals.

        Example:
            >>> canonical = mol.electronic_state.canonical
            >>> atomic = mol.electronic_state.basis
            >>> overlaps = canonical.overlap(atomic)
            >>> overlaps[i, j] == canonical.orbitals[i].overlap(atomic.orbitals[j])
            True
        """
        return self.dot(self.aobasis.overlaps.dot(other.T))

    def __iter__(self):
        return iter(self.orbitals)

    def __len__(self):
        return len(self.orbitals)

    def __str__(self):
        return '%s orbitals' % self.orbtype

    def __repr__(self):
        return '<%d %s %s in %s>' % (self.num_orbs, self.orbtype,
                                     self.__class__.__name__, str(self.wfn))

    @property
    def energies(self):
        """u.Vector[energy]: energies of the molecular orbitals

        This is just the diagonal of the fock matrix"""
        return self.fock.diagonal()

    @property
    def fock(self):
        """u.Array[energy]: Fock matrix for these orbitals"""
        return self.from_ao(self.wfn.fock_ao)

    @property
    def overlaps(self):
        """np.array: overlap matrix for these orbitals"""
        return self.from_ao(self.wfn.aobasis.overlaps)

    @property
    def h1e(self):
        """u.Array[energy]: 1-electron matrix elements for these orbitals"""
        return self.from_ao(self.wfn.aobasis.h1e)

    @property
    def h2e(self):
        """u.Array[energy]: 2-electron matrix elements for these orbitals"""
        return self.fock - self.h1e

    def from_ao(self, ao_operator):
        """ Transform an operator into this orbital basis from the ao basis

        Given the matrix elements :math:`\hat O_{\mu \nu}` of an operator over AO basis indices
        :math:`\mu,\nu`, returns the operator's matrix elements :math:`\hat O_{ij}` over
        orbital indices :math:`i,j`:

        ..math::
            \hat O_{ij} =
             \left \langle i \right| \hat O \left| j \right \rangle =
              \sum_{\mu \nu}C_{i \mu} O_{\mu \nu} C_{j \nu}


        where :math:`C_{i \mu}` is the expansion coefficient for AO basis function :math:`\mu` in
        molecular orbital _i_.

        Args:
            ao_operator (u.Array): matrix elements of the operator in the ao basis

        Returns:
            u.Array: matrix elements of the operator in this orbital basis

        Note:
            Assumes that this set of orbitals is orthogonal
        """
        # Dot doesn't place nice with units, so we need to pass them explicitly
        ao_units = u.get_units(ao_operator)
        return self.coeffs.dot(ao_operator.dot(self.coeffs.T)) * ao_units

    def to_ao(self, mo_operator):
        """ Transform an operator from this orbital basis into the AO basis

        Given the matrix elements :math:`\hat O_{ij}` of an operator over orbital basis indices
        :math:`i,j`, returns the operator's matrix elements :math:`\hat O_{\mu \nu}` over
        orbital indices :math:`\mu, \nu`:

        ..math::
           \hat O_{\mu \nu} =
           \left \langle \mu \right| \hat O \left| \nu \right \rangle =
           \sum_{i,j,\lambda,\kappa}S_{\mu \lambda} C_{i \lambda} O_{ij} C_{j \kappa} S_{\kappa \nu}

        where :math:`S_{\mu \nu} = \left \langle \mu | \nu \right \rangle` is the AO overlap matrix
        and :math:`C_{i \mu}` is the expansion coefficient for AO basis function :math:`\mu` in
        molecular orbital _i_.

        Args:
            mo_operator (u.Array): matrix elements of the operator in this orbital basis

        Returns:
            u.Array: matrix elements of the operator in the AO basis
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


class Attribute(object):
    """For overriding a property in a superclass - turns the attribute back
    into a normal instance attribute"""
    def __init__(self, name):
        self.name = name

    def __get__(self, instance, cls):
        return getattr(instance, self.name)

    def __set__(self, instance, value):
        return setattr(instance, self.name, value)


class BasisSet(MolecularOrbitals):
    """
    Stores a basis, typically of atomic orbitals.

    This is a special orbital type
    """
    overlaps = Attribute('_overlaps')
    h1e = Attribute('_h1e')

    def __init__(self, mol, orbitals, name=None,
                 h1e=None, overlaps=None,
                 angulartype=None,
                 **kwargs):
        self.mol = mol
        self.orbitals = orbitals
        self.coeffs = np.identity(len(self.orbitals))
        self.basis = self
        self.orbtype = 'aobasis'
        self.angulartype = angulartype
        assert self.angulartype in (None, 'spherical', 'cartesian')

        self.wfn = None
        for iorb, orbital in enumerate(self.orbitals):
            orbital.index = iorb
            orbital.coeffs = self.coeffs[iorb, :]

        self.basisname = name
        self.h1e = h1e
        self.overlaps = overlaps
        for kw, val in kwargs.iteritems():
            setattr(self, kw, val)

        self.on_atom = {}
        for fn in self.orbitals:
            self.on_atom.setdefault(fn.atom, []).append(fn)

    def calculate_orb_grid(self, coeffs, padding=2.5*u.angstrom, npoints=30):
        """Abstract method - subclasses need to implement this"""
        raise NotImplementedError
        # TODO - generic implementation here
        # TODO: needs to go fast for gaussians - use numpy vectorization
        # TODO: balance caching, precomputing with mem usage ...

    def __repr__(self):
        return '<%s (%s) of %s>' % (self.__class__.__name__, self.basisname, self.mol)


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
            self.orbitals['atomic'] = self.aobasis

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
