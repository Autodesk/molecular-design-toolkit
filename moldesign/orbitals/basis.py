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


@toplevel
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
            orbital.basis = self

        self.basisname = name
        self.h1e = h1e
        self.overlaps = overlaps
        for kw, val in kwargs.items():
            setattr(self, kw, val)

        self.on_atom = {}
        for fn in self.orbitals:
            self.on_atom.setdefault(fn.atom, []).append(fn)

    def __call__(self, coords, coeffs=None):
        """Calculate the value of each basis function at the specified coordinates.

        Note:
            This is just a pass-through to a specific implementation - PYSCF's eval_ao routine
            for now.

        Args:
            coords (Array[length]): List of coordinates.
            coeffs (Vector): List of ao coefficients (optional; if not passed, all basis fn
                 values are returned)

        Returns:
            Array[length]: if ``coeffs`` is not passed, an array of basis fn values at each
               coordinate. Otherwise, a list of orbital values at each coordinate
        """
        from moldesign.interfaces.pyscf_interface import basis_values
        return basis_values(self.wfn.mol, self, coords, coeffs=coeffs,
                            positions=self.wfn.positions)

    def __repr__(self):
        return '<%s (%s) of %s>' % (self.__class__.__name__, self.basisname, self.mol)

    @property
    def fock(self):
        return self.wfn.fock_ao

    @property
    def density_matrix(self):
        return self.wfn.density_matrix_ao

