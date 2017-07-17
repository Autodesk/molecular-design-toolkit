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

from ..utils import Attribute
from . import toplevel, MolecularOrbitals


@toplevel
class BasisSet(MolecularOrbitals):
    """
    Stores a basis, typically of atomic orbitals.

    This is a case of a general set of molecular orbitals, where the coefficient matrix
    must be the identity

    Args:
        mol (mdt.Molecule): the molecule these basis functions belong to
        orbitals (List[AtomicBasisFunction]): list of basis functions comprising this set
        name (str): name of this basis set
        h1e (Matrix[energy]): 1-electron elements of the hamiltonian
        overlaps (Matrix): overlap matrix
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
            orbital.basis = self

        self.basisname = name
        self.h1e = h1e
        self.overlaps = overlaps
        for kw, val in kwargs.items():
            setattr(self, kw, val)

        self._basis_fn_by_atom = {}
        for fn in self.orbitals:
            self._basis_fn_by_atom.setdefault(fn.atom.index, []).append(fn)

    def __call__(self, coords, coeffs=None):
        """Calculate the amplitude of a list of orbitals at the specified coordinates.

        Returns an array of orbital amplitudes. The amplitude of the _n_th orbital at the
        _j_th coordinate is stored at position ``(j, n)`` in the returned array.

        Args:
            coords (Matrix[shape=(*,3)]): List of coordinates.
            coeffs (Matrix[shape=(*, nbasis)]): List of ao coefficients (optional;
               if not passed, the amplitudes of all basis functions will be returned

        Returns:
            Matrix:
               - if ``coeffs`` is passed, an array of orbital amplitudes at the given positions
                  of size ``(len(coords), len(coeffs))``.
               - if ``coeffs`` is NOT passed, an array of basis function amplitudes
                 of size ``(len(coords), len(aobasis))``.
        """
        basis_vals = np.zeros((len(coords), len(self))) * self.orbitals[0]._get_wfn_units()
        for ibf, bf in enumerate(self.orbitals):
            basis_vals[:, ibf] = bf(coords)

        if coeffs is None:
            return basis_vals
        else:
            return np.dot(basis_vals, coeffs.T)

    def get_basis_functions_on_atom(self, atom):
        """ Return a list of basis functions on this atom

        Args:
            atom (moldesign.Atom): query atom

        Returns:
            List[AtomicBasisFunction]: basis functions centered on this atom
        """
        return self._basis_fn_by_atom[atom.index]

    def __repr__(self):
        return '<%s (%s) of %s>' % (self.__class__.__name__, self.basisname, self.mol)

    @property
    def fock(self):
        return self.wfn.fock_ao

    @property
    def density_matrix(self):
        return self.wfn.density_matrix_ao

