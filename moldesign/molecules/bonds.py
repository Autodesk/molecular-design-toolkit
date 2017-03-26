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
        if a2 in self._graph[a1]:
            raise ValueError("%s and %s already bonded, can't create another bond" %
                             (a1, a2))
        bond = Bond(a1, a2)
        bond.order = order
        return bond

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
        for nbr, order in atom.bond_graph.iteritems():
            if nbr.molecule is self._mol:
                self._graph[atom][nbr] = self._graph[nbr][atom] = order

    def _removeatom(self, atom):
        """ Private data mangement method - removes this atom from the bond graph.
        """
        for nbr in self._graph[atom]:
            self._graph[nbr].pop(atom)
        self._graph.pop(atom)

    def __getitem__(self, atoms):
        if atoms[0] not in self._graph or atoms[1] not in self._graph[atoms[0]]:
            raise ValueError('%s is not bonded to %s' % (atoms[0], atoms[1]))
        return Bond(atoms[0], atoms[1])


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

        if not a1.bonded_to(a2):
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
