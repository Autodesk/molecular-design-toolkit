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

from . import toplevel


@toplevel
class Bond(object):
    """
    A bond between two atoms.

    Args:
        a1 (Atom): First atom
        a2 (Atom): Second atom (the order of atoms doesn't matter)
        order (int): Order of the bond

    Notes:
        Comparisons and hashes involving bonds will return True if the atoms involved in the bonds
        are the same. Bond orders are not compared.

        These objects are used to represent and pass bond data only - they are not used for storage.

    Attributes:
        a1 (Atom): First atom in the bond; assigned so that ``self.a1.index < self.a2.index``
        a2 (Atom): Second atom in the bond; assigned so that ``self.a2.index > self.a1.index``
        order (int): bond order (can be ``None``); not used in comparisons
    """
    def __init__(self, a1, a2, order=None):
        if a1.molecule is not a2.molecule:
            raise ValueError('Cannot create bond for atoms in different molecules.')
        else:
            self.molecule = a1.molecule

        if a1.index is not None and a2.index is not None and a1.index > a2.index:
            a1, a2 = a2, a1
        self.a1 = a1
        self.a2 = a2
        if order is None:
            try: self.order = self.a1.bond_graph[a2]
            except KeyError: self.order = None
        else:
            self.order = order

    def __eq__(self, other):
        return (self.a1 is other.a1) and (self.a2 is other.a2)

    def __hash__(self):
        """Has this object using the atoms involved in its bond"""
        return hash((self.a1, self.a2))

    def partner(self, atom):
        """ Return this atom's *partner* in the bond -- i.e., the other atom in the bond

        Args:
            atom (mdt.Atom): return the atom that this one is bonded to

        Returns:
            mdt.Atom: the passed atom's partner

        Raises:
            ValueError: if the passed atom is not part of this bond
        """
        if atom is self.a1:
            return self.a2
        elif atom is self.a2:
            return self.a1
        else:
            raise ValueError('%s is not part of this bond' % atom)

    @property
    def length(self):
        return self.a1.distance(self.a2)

    @property
    def name(self):
        """ str: name of the bond """
        return '{a1.name} (#{a1.index}) - {a2.name} (#{a2.index}) (order: {order})'.format(
                        a1=self.a1, a2=self.a2, order=self.order)

    @property
    def ff(self):
        """mdt.forcefield.BondTerm: the force-field term for this bond (or ``None`` if no
            forcefield is present)
        """
        return self.molecule.ff.get_bond_term(self)
