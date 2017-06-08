""" Class definitions for atomic and molecular orbitals.

Notes:
    In this documentation, we use the following conventions for labeling orbitals:
      - atomic orbitals using lower case greek labels and subscripts, e.g.,
          :math:`\left| \mu \right \rangle, F_{\nu \lambda}, etc.
      - molecular orbitals use lower case labels and subscripts, e.g.,
          :math:`\left| i \right \rangle, F_{kl}, etc.
      - adiabatic electronic states are indexed using capital letters, _N_, _L_, _M_, etc.

"""

from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library
standard_library.install_aliases()

# Copyright 2017 Autodesk Inc.
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

import numpy as np

from moldesign import units as u
from moldesign.utils import Alias

SHELLS = {0: 's', 1: 'p', 2: 'd', 3: 'f', 4: 'g', 5: 'h'}
ANGMOM = {v: k for k, v in SHELLS.items()}

# See https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
SPHERICALNAMES = {(0, 0): 's', (1, -1): 'p(x)', (1, 0): 'p(z)', (1, 1): 'p(y)',
                  (2, 0): 'd(z^2)', (2, -2): 'd(xy)', (2, -1): 'd(yz)',
                  (2, 1): 'd(xz)', (2, 2): 'd(x^2-y^2)',
                  (3, -3): 'f(3yx^2-y^3)', (3, -2): 'f(xyz)', (3, -1): 'f(yz^2)',
                  (3, 0): 'f(z^3)', (3, 1): 'f(xz^2)', (3, 2): 'f(zx^2-zy^2)',
                  (3, 3): 'f(x^3-3xy^2)'}
ANGULAR_NAME_TO_COMPONENT = {'': (0, 0), 'x': (1, -1), 'y': (1, 1), 'z': (1, 0),
                             'z^2': (2, 0), 'xy': (2, -2), 'yz': (2, -1), 'xz': (2, 1),
                             'x^2-y^2': (2, 2),
                             'zx^2-zy^2': (3, 2), 'xyz': (3, -2), 'z^3': (3, 0),
                             '3yx^2-y^3': (3, -3), 'x^3 - 3xy^2': (3, 3), 'xz^2': (3, 1),
                             'yz^2': (3, -1)}


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
                 occupation=None, name='unnamed'):
        """ Initialization:

        Args:
            coeffs (numpy.array): orbital coefficients
            basis (moldesign.orbitals.basis.BasisSet): basis for this orbital
            wfn (moldesign.orbitals.wfn.ElectronicWfn): total electronic wavefunction that this orbital is a part of
            occupation (float): occupation of this orbital
            name (str): optional label for this orbital:
        """
        self.coeffs = np.array(coeffs)
        self.name = name
        self.molecule = None
        self.basis = basis
        self.occupation = occupation
        self.wfn = wfn
        self.index = None  # will be set by the containing MolecularOrbitals object

        # Assign the basis functions
        if wfn is not None and wfn.aobasis is not None:
            if self.basis is not None:
                assert wfn.aobasis is self.basis
            else:
                self.basis = wfn.aobasis
        if self.basis is not None:
            assert len(self.coeffs) == len(self.basis)

    def overlap(self, other):
        """ Calculate overlap with another orbital

        Args:
            other (Orbital): calculate the overlap with this orbital

        Returns:
            float: orbital overlap
        """
        return self.coeffs.dot(self.basis.overlaps.dot(other.coeffs))

    def fock_element(self, other):
        """ Calculate fock matrix element with another orbital

        Args:
            other (Orbital): calculate the fock element with this orbital

        Returns:
            u.Scalar[energy]: fock matrix element
        """
        return self.wfn.fock_ao.dot(other.coeffs).ldot(self.coeffs)  # uses ldot to preserve units

    @property
    def energy(self):
        """ u.Scalar[energy]: This orbital's energy

        Note:
            This is equivalent to self.fock(self)
        """
        return self.fock_element(self)

    def __call__(self, coords):
        """ Calculate the orbital's value at a position (or list of positions)

        Args:
            coords (u.Vector[length]): Coordinates (shape ``(len(coords), 3)``)

        Returns:
            u.Scalar[length**(-3/2)]: value of the orbital
        """
        return self.basis(coords, coeffs=self.coeffs)

    def __repr__(self):
        return '<%s %s>' % (self.__class__.__name__, self.name)

    __len__ = Alias('coeffs.__len__')


class MolecularOrbitals(object):
    """
    Stores a wfn of molecular orbitals in an AO wfn
    Orbitals are accessed as orbs[orbital index, ao index]
    """
    def __init__(self, orbitals, wfn=None, basis=None,
                 canonical=False,
                 orbtype=None):
        """ Initialization:

        Args:
            orbitals (List[Orbital] OR numpy.array): EITHER a list of orbitals OR an
               array of coefficients (indexed as coeffs[orbital_index, basis_index])
            basis (moldesign.orbitals.basis.BasisSet): orbital basis set
            wfn (moldesign.orbitals.wfn.ElectronicWfn): total electronic wavefunction that this orbital is a part of
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
            if thisorb.overlap(otherorb) < -1.0 * threshold:
                thisorb.coeffs *= -1.0
                # TODO: print a warning if overlap is small?

    def overlap(self, other):
        """ Calculate overlaps between this and another set of orbitals

        Args:
            other (MolecularOrbitals):

        Returns:
            numpy.ndarray: overlaps between the two sets of orbitals

        Example:
            >>> canonical = mol.wfn.canonical
            >>> atomic = mol.wfn.basis
            >>> overlaps = canonical.overlap(atomic)
            >>> overlaps[i, j] == canonical.orbitals[i].overlap(atomic.orbitals[j])
            True
        """
        return self.coeffs.dot(self.wfn.aobasis.overlaps.dot(other.coeffs.T))

    def __iter__(self):
        return iter(self.orbitals)

    def __len__(self):
        return len(self.orbitals)

    def __str__(self):
        return '%s orbitals' % self.orbtype

    def __repr__(self):
        return '<%d %s %s in %s>' % (len(self), self.orbtype,
                                     self.__class__.__name__, str(self.wfn))

    def __getitem__(self, item):
        return self.orbitals[item]

    def _to_ao_density_matrix(self):
        c = self.coeffs * self.occupations[:, None]/2.0
        return 2.0*c.T.dot(c)

    @property
    def energies(self):
        """u.Vector[energy]: energies of the molecular orbitals

        This is just the diagonal of the fock matrix"""
        return self.fock.diagonal()

    @property
    def occupations(self):
        """ np.ndarray: orbital occupation numbers
        """
        return np.array([orb.occupation for orb in self.orbitals])

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


