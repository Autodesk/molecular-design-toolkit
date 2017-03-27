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

import moldesign
from moldesign import data, utils


class BondGraph(object):
    def __init__(self, mol):
        self._mol = mol
        self._graph = {}

    def __iter__(self):
        for atom in self._graph:
            for nbr in self._graph[atom]:
                if atom.index > nbr.index:
                    continue  # don't double count
                yield Bond(atom, nbr)

    def create(self, a1, a2, order):
        """ Create a bond between two atoms in this molecule

        Args:
            a1 (Atom): first atom in bond
            a2 (Atom): second atom in bond
            order (int): bond order

        Returns:
            Bond: record of the bond
        """
        if a1 is a2:
            raise ValueError("Cannot bond atom to itself!")

        if a2 in self._graph[a1]:
            raise ValueError("%s and %s already bonded, can't create another bond" %
                             (a1, a2))

        self._graph[a1][a2] = order
        self._graph[a2][a1] = order
        return Bond(a1, a2)

    def delete(self, atom_or_bond, a2=None):
        """ Delete a bond

        Args:
            atom_or_bond (Atom or Bond): first atom in bond (or bond to delete
            a2 (Atom): second atom in bond

        Returns:
            Bond: record of the bond (which is no longer part of the molecule)
        """
        if a2 is None:
            b = atom_or_bond
        else:
            b = Bond(atom_or_bond, a2)
        if b not in self:
            raise ValueError("Bond %s is not part of this molecule")
        b.delete()
        return b

    def __contains__(self, bond):
        return bond.a1 in self._graph and bond.a2 in self._graph[bond.a1]

    def _addatom(self, atom):
        """ Private data mangement method - adds entry for this atom in the bond graph.
        Only adds bonds for atoms that are already tracked in this graph.
        """
        assert atom not in self._graph
        self._graph[atom] = {}
        for nbr, order in atom._graph.iteritems():
            if nbr.molecule is self._mol:
                self._graph[atom][nbr] = order
                if nbr in self._graph:
                    self._graph[nbr][atom] = order

    def _removeatom(self, atom):
        """ Private data mangement method - removes this atom from the bond graph.
        """
        for nbr in self._graph[atom]:
            self._graph[nbr].pop(atom)
        self._graph.pop(atom)

    def __getitem__(self, atoms):
        try:
            a1, a2 = atoms
        except TypeError:
            return AtomBonds(atoms, self._graph[atoms])
        else:
            return Bond(a1, a2)


class AtomBonds(object):
    """ A view of one atom's bonds
    """
    __len__ = utils.Alias('_graph.__len__')

    def __init__(self, atom, graphdict):
        self.atom = atom
        self._graph = graphdict

    def __getitem__(self, atom):
        return Bond(self.atom, atom)

    def __iter__(self):
        for nbr in self._graph:
            yield Bond(self.atom, nbr)

    def __contains__(self, other):
        if isinstance(other, moldesign.Atom):
            return other in self._graph
        elif isinstance(other, Bond):
            if other.a1 is self:
                return other.a2 in self
            elif other.a2 is self:
                return other.a1 in self
            else:
                return False
        else:
            return False

    def iteritems(self):
        """ Yield tuples of this atom's bonded neighbors and the corresponding bond orders

        Yields:
            Tuple[moldesign.Atom, int]: ([atom bonded to this one], [order of the bond])
        """
        return self._graph.iteritems()

    @property
    def heavy(self):
        """ Iterator[Bond]: list of all heavy atom bonds (where BOTH atoms are not hydrogen)

        Note:
            this returns an empty list if called on a hydrogen atom
        """
        for bond in self:
            if bond.a1.atnum > 1 and bond.a2.atnum > 1:
                yield bond

    @property
    def atoms(self):
        """ List[Atom]: list of atoms this atom is bonded to
        """
        return self._graph.keys()

    def create(self, other, order):
        """ Create or modify a bond with another atom

        Args:
            other (Atom): atom to bond to
            order (int): bond order

        Returns:
            moldesign.molecules.bonds.Bond: bond object
        """
        if other is self:
            raise ValueError("Cannot bond atom to itself!")
        if self.atom.molecule and (self.atom.molecule is other.molecule):
            self.atom.molecule.bonds.create(self.atom, other, order)
        else:  # allow unassigned atoms to be bonded to anything for building purposes
            if self.atom.molecule is None:
                self._graph[other] = order
            if other.molecule is None:
                other._graph[self.atom] = order
            if self.atom.molecule is not None and other.molecule is not None:
                raise ValueError("Can't bond atom in different molecules")
        return Bond(self.atom, other)

    def delete(self, other):
        """ Delete the bond between this atom and another

        Args:
            other (moldesign.Atom): Atom to delete the bond with

        Raises:
            ValueError: if the two atoms are not bonded
        """
        if self.atom.molecule is other.molecule:
            self.atom.molecule.bonds.delete(self, other)
        else:
            if self.atom.molecule is None:
                self._graph.pop(other)
            if other.molecule is None:
                other._graph.pop(self, None)
            if self.atom.molecule is not None and other.molecule is not None:
                raise ValueError("Can't bond atom in different molecules")

    def __str__(self):
        start = "%s bonds to:" % self.atom
        s = ', '.join('%s (%s)' % (nbr, data.BONDNAMES.get(order, 0))
                      for nbr, order in self._graph.iteritems())
        return start + s

    def __repr__(self):
        return '<%s>' % self


class Bond(object):
    """
    A bond between two atoms.

    Args:
        a1 (Atom): First atom
        a2 (Atom): Second atom (the order of atoms doesn't matter)

    Notes:
        Comparisons and hashes involving bonds will return True if the atoms involved in the bonds
        are the same. Bond orders are not compared.

        These objects are used to represent and pass bond data only - they are not used for storage.

    Raises:
        ValueError: if these two atoms are not bonded

    Attributes:
        a1 (Atom): First atom in the bond; assigned so that ``self.a1.index < self.a2.index``
        a2 (Atom): Second atom in the bond; assigned so that ``self.a2.index > self.a1.index``
        order (int): bond order (can be ``None``); not used in comparisons
    """
    def __init__(self, a1, a2):
        if a1.molecule is not a2.molecule:
            raise ValueError('Cannot create bond for atoms in different molecules.')
        else:
            self.molecule = a1.molecule

        if a2 not in a1.bonds:
            raise ValueError('%s is not bonded to %s' % (a1, a2))

        self._exists = True

        if a1.index > a2.index:
            a1, a2 = a2, a1
        self.a1 = a1
        self.a2 = a2

    def __eq__(self, other):
        return (self.a1 is other.a1) and (self.a2 is other.a2)

    @property
    def order(self):
        if not self._exists:
            raise ValueError("This bond has been deleted.")
        return self.molecule.bonds._graph[self.a1][self.a2]

    @order.setter
    def order(self, o):
        self.molecule.bonds._graph[self.a1][self.a2] = o
        self.molecule.bonds._graph[self.a2][self.a1] = o
        self._exists = True

    def delete(self):
        self.molecule.bonds._graph[self.a1].pop(self.a2)
        self.molecule.bonds_graph[self.a2].pop(self.a1)
        self._exists = False

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
        if self._exists:
            return '{a1.name} (#{a1.index}) - {a2.name} (#{a2.index}) (order: {order})'.format(
                            a1=self.a1, a2=self.a2, order=self.order)
        else:
            return '{a1.name} (#{a1.index}) - {a2.name} (#{a2.index}) (deleted)'.format(
                            a1=self.a1, a2=self.a2)

    @property
    def ff(self):
        """mdt.forcefield.BondTerm: the force-field term for this bond (or ``None`` if no
            forcefield is present)
        """
        return self.molecule.ff.get_bond_term(self)
