# coding: utf-8

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
from past.builtins import basestring

import numpy as np

from . import toplevel
from .. import units as u
from .. import mathutils


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
    SYMBOLS = {1: u'-', 2: u'=', 3: u'≡'}

    def __init__(self, a1, a2, order=None):
        if a1.molecule is not a2.molecule:
            raise ValueError('Cannot create bond for atoms in different molecules.')

        if a1.index is not None and a2.index is not None and a1.index > a2.index:
            a1, a2 = a2, a1
        self.a1 = a1
        self.a2 = a2
        if order is None:
            try: self.order = self.a1.bond_graph[a2]
            except KeyError: self.order = None
        else:
            self.order = order

    def __str__(self):
        return "%s bond between %s and %s (order: %s)" % (self.type,
                                                          self.a1._shortstr(), self.a2._shortstr(),
                                                          self.order)

    def __repr__(self):
        return "<%s>" % self

    @property
    def type(self):
        return self.a1.symbol + self.SYMBOLS.get(self.order, u'?̶') + self.a2.symbol

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
    def molecule(self):
        """moldesign.Molecule: molecule this bond belongs to (or None if not assigned

        Raises:
            ValueError: if the atoms are assigned to different molecules
        """
        if self.a1.molecule is not self.a2.molecule:
            raise ValueError("Atoms in %s belong to different molecules" % self)
        else:
            return self.a1.molecule

    @property
    def midpoint(self):
        return (self.a1.position + self.a2.position) / 2.0

    def align(self, other, centered=True):
        """ Rotates the entire molecule to align this bond with another object.

        Args:
            other (str or Vector[len=3] or Bond): Object to align this bond with, which may be:
               - a string: 'x', 'y', or 'z',
               - a len-3 vector, or
               - another :class:`Bond` object
            centered (bool): if True (default), center this bond at the origin OR the midpoint of
             the other bond
        """
        mol = self.molecule

        centering = -self.midpoint

        if isinstance(other, Bond):
            direction = (other.a2.position - other.a1.position).normalized()
            if centered:
                centering += other.midpoint

        elif isinstance(other, basestring):
            arr = np.zeros(3)
            arr[DIMLABELS[other]] = 1.0
            direction = arr

        else:
            direction = other

        target = mathutils.normalized(u.array(direction))
        vec = (self.a2.position - self.a1.position).normalized()

        angle, normal = mathutils.alignment_rotation(vec, target)

        if centered:
            mol.positions += centering

        if abs(mathutils.normalize_angle(angle)) > 1e-3 * u.degrees:
            mol.rotate(angle, normal, center=self.midpoint)

    @property
    def ff(self):
        """mdt.forcefield.BondTerm: the force-field term for this bond (or ``None`` if no
            forcefield is present)
        """
        if self.molecule.ff is not None:
            return self.molecule.ff.get_bond_term(self)
        else:
            return None

DIMLABELS = {'x': 0, 'y': 1, 'z':2}
