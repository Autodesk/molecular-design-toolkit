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
from . import PrimitiveSum, SHELLS, SPHERICALNAMES


class AtomicBasisFunction(PrimitiveSum):
    """ Stores an atomic basis function.

    Note:
        Either l and m should be passed, or cart, but not both.

    Args:
        atom (moldesign.Atom): The atom this basis function belongs to
        index (int): the index of this basis function (it is stored as
            ``wfn.basis[self.index]``)
        n (int): principal quantum number (``n>=1``) - shell metadata (used for labeling only)
        l (int): total angular momentum quantum number (``l<=n-1``)
        m (int): z-angular momentum quantum number (optional -
             for spherical sets only; ``|m|<=l``)
        cart (str): cartesian component (optional; for cartesian sets only)
        primitives (List[PrimitiveBase]): List of primitives, if available
    """
    def __init__(self, atom, n=None, l=None, m=None, cart=None, primitives=None):
        super().__init__(primitives)
        self._atom = atom
        self.atom_index = atom.index
        self.n = n
        self.l = l
        self.m = m
        if cart is not None:
            assert self.m is None, 'Both cartesian and spherical components passed!'
            assert len(cart) == self.l, \
                'Angular momentum does not match specified component %s' % cart
            for e in cart:
                assert e in 'xyz'
            self.cart = ''.join(sorted(cart))
        else:
            self.cart = None

        # These quantities can't be defined until we assemble the entire basis
        self.coeffs = None
        self.molecule = atom.molecule
        self.basis = None
        self.occupation = None
        self.wfn = None

    @property
    def atom(self):
        """ moldesign.Atom: the atom this basis function belongs to

        We get the atom via an indirect reference, making it easier to copy the wavefunction
        """
        if self.wfn is not None:
            return self.wfn.mol.atoms[self.atom_index]
        else:
            return self._atom

    @property
    def orbtype(self):
        """ A string describing the orbital's angular momentum state.

        Examples:
            >>> AtomicBasisFunction(n=1, l=0).orbtype
            's'
            >>> AtomicBasisFunction(n=2, l=1, cart='y').orbtype
            'py'
            >>> AtomicBasisFunction(n=3, l=2, m=0).orbtype
            'd(z^2)'
        """
        if self.l == 0: t = 's'
        elif self.cart is not None: t = SHELLS[self.l] + self.cart
        else: t = SPHERICALNAMES[self.l, self.m]
        return t

    @property
    def aotype(self):
        """ A string describing the orbital's state.

        Examples:
            >>> AtomicBasisFunction(n=1, l=0).aotype
            '1s'
            >>> AtomicBasisFunction(n=2, l=1, cart='y').aotype
            '2py'
            >>> AtomicBasisFunction(n=3, l=2, m=0).aotype
            '3d(z^2)'
        """
        t = self.orbtype
        if self.n:
            return '%s%s' % (self.n, t)
        else:
            return t

    def __str__(self):
        return 'AO ' + self.name

    @property
    def name(self):
        try:
            return '%s on atom %s' % (self.aotype, self.atom.name)
        except:
            return 'Basis Fn'

    def __repr__(self):
        return '<%s %s>' % (self.__class__.__name__, self.name)