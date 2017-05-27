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
import copy
import collections
import itertools

import numpy as np

import moldesign as mdt
from .. import units as u
from .. import utils, external, mathutils, helpers
from . import toplevel


class AtomGroup(object):
    """ Mixin functions for objects that have a ``self.atoms`` attribute with a list of atoms

    Attributes:
        atoms (List[Atom]): a list of atoms
    """
    draw2d = helpers.WidgetMethod('atomgroups.draw2d')
    draw3d = helpers.WidgetMethod('atomgroups.draw3d')
    draw = helpers.WidgetMethod('atomgroups.draw')

    def __init__(self, *args, **kwargs):
        """ This should never be called directly - it will be called by the `super` methods
        of its subclasses """
        super().__init__(*args, **kwargs)
        self._atom_attrs = None
        self.viz2d = None
        self.viz3d = None

    @property
    def num_atoms(self):
        """ int: number of atoms in this object """
        return len(self.atoms)
    natoms = num_atoms

    @property
    def heavy_atoms(self):
        """ AtomList: a list of all heavy atoms (i.e., non-hydrogen) in this object """
        return AtomList([a for a in self.atoms if a.atnum != 1])

    @property
    def mass(self):
        """ u.Scalar[mass]: total mass of this object
        """
        return sum(a.mass for a in self.atoms)

    @property
    def momentum(self):
        """ u.Vector[momentum]: total momentum of this object
        """
        return self.momenta.sum(axis=0)

    @property
    def velocity(self):
        """ u.Vector[velocity]: center of mass velocity of this object
        """
        return old_div(self.momentum,self.mass)

    @property
    def kinetic_energy(self):
        r""" u.Scalar[energy]: Classical kinetic energy :math:`\sum_{\text{atoms}} \frac{p^2}{2m}`
        """
        return helpers.kinetic_energy(self.momenta, self.masses)

    def get_atoms(self, *keywords, **queries):
        """Allows keyword-based atom queries. Returns atoms that match ALL queries.

        Args:
            *keywords (list): pre-set keywords (currently, just selects by residue type)
            **queries (dict): attributes (or residue attributes) to match

        Examples:
            >>> mol.get_atoms('protein')  # returns all atoms in proteins
            >>> mol.get_atoms(name='CA')  # returns all alpha carbons
            >>> mol.get_atoms('dna', symbol='P')  # returns all phosphorus in DNA

        Returns:
            AtomList: the atoms matching this query
        """
        if not (queries or keywords):
            return mdt.AtomList(self.atoms)

        atoms = self.atoms

        for key in keywords:
            if key in 'protein dna rna water unknown':
                atoms = mdt.AtomList(atom for atom in atoms
                                     if atom.residue.type == key)

        result = mdt.AtomList()
        for atom in atoms:
            for field, val in queries.items():
                if getattr(atom, field, None) != val and getattr(atom.residue, field, None) != val:
                    break
            else:
                result.append(atom)

        return result

    def calc_distance_array(self, other=None):
        """ Calculate an array of pairwise distance between all atoms in self and other

        Args:
            other (AtomContainer): object to calculate distances to (default: self)

        Returns:
            u.Array[length]: 2D array of pairwise distances between the two objects

        Example:
            >>> dists = self.calc_distance_array(other)
            >>> dists[i, j] == self.atoms[i].distance(other.atoms[j])
        """
        from scipy.spatial.distance import cdist

        other = utils.if_not_none(other, self)
        try:
            other_positions = other.positions.defunits_value()
        except AttributeError:
            other_positions = np.array([other.position.defunits_value()])

        distances = cdist(self.positions.defunits_value(), other_positions)
        return distances * u.default.length

    def calc_displacements(self):
        """ Calculate an array of displacements between all atoms in this object

        Returns:
            u.Array[length]: array of pairwise displacements between atoms

        Example:
            >>> displacements = self.calc_displacements(other)
            >>> displacements[i, j] == (self.atoms[i].position - self.atoms[j].position)
        """
        # TODO: allow other, similar to calc_distance array
        return utils.pairwise_displacements(self.positions)

    def distance(self, other):
        """Returns closest distance between this and the other entity

        Args:
            other (AtomContainer): object to calculate distance to

        Returns:
            u.Scalar[length]: closest distance between self and other

        Example:
            >>> distance = self.distance(other)
            >>> distance == self.calc_distance_array(other).min()
        """
        distance_array = self.calc_distance_array(other)
        return distance_array.min()

    @property
    def center_of_mass(self):
        """ units.Vector[length]: The (x,y,z) coordinates of this object's center of mass """
        if self.num_atoms == 0:  # nicer exception than divide-by-zero
            raise ValueError('"%s" has no atoms' % str(self))

        total_mass = 0.0 * u.default.mass
        com = np.zeros(3) * u.default.length * u.default.mass
        for atom in self.atoms:
            total_mass += atom.mass
            com += atom.position * atom.mass
        com = com / total_mass
        return com

    @center_of_mass.setter
    def center_of_mass(self, value):
        vec = value - self.com
        self.translate(vec)

    com = center_of_mass  # synonym

    def _getatom(self, a):
        """ Given an atom's name, index, or object, return the atom object
        """
        if a is None:
            return None
        elif isinstance(a, basestring) or isinstance(a, int):
            return self[a]
        else:
            return a

    def angle(self, a1, a2, a3):
        """ Calculate the angle between three atoms.

        Atoms can be passed as the atoms themselves or as the atom names

        Args:
            a1, a2, a3 (str OR int OR moldesign.Atom): atoms defining the angle

        Returns:
            units.Scalar[angle]
        """
        # TODO: use single dispatch to also accept two bonds
        return mdt.geom.angle(*list(map(self._getatom, (a1, a2, a3))))

    def dihedral(self, a1, a2, a3=None, a4=None):
        """ Calculate the dihedral angle between atoms a1, a2, a3, a4.

        Atoms can be passed as the atoms themselves or as the atom names

        Args:
            a1, a2, a3, a4 (str OR int OR moldesign.Atom): atoms defining the dihedral

        Returns:
            units.Scalar[angle]
        """
        return mdt.geom.dihedral(*list(map(self._getatom, (a1, a2, a3, a4))))

    def copy_atoms(self):
        """
        Copy a group of atoms which may already have bonds, residues, and a parent molecule
        assigned. Do so by copying only the relevant entities, and creating a "mask" with
        deepcopy's memo function to stop anything else from being copied.

        Returns:
            AtomList: list of copied atoms
        """
        from . import ChildList

        oldatoms = self.atoms
        old_bond_graph = {a: {} for a in self.atoms}
        for atom in self.atoms:
            for nbr in atom.bond_graph:
                if nbr in old_bond_graph:
                    old_bond_graph[atom][nbr] = atom.bond_graph[nbr]

        # Figure out which bonds, residues and chains to copy
        tempatoms = AtomList([copy.copy(atom) for atom in oldatoms])
        old_to_new = dict(zip(oldatoms, tempatoms))
        temp_bond_graph = {a: {} for a in tempatoms}
        replaced = {}
        for atom, oldatom in zip(tempatoms, oldatoms):
            atom.molecule = None
            atom.bond_graph = {}

            if atom.chain is not None:
                if atom.chain not in replaced:
                    chain = copy.copy(atom.chain)
                    chain.molecule = None
                    chain.children = ChildList(chain)
                    replaced[atom.chain] = chain
                else:
                    chain = replaced[atom.chain]

                atom.chain = chain

            if atom.residue is not None:
                if atom.residue not in replaced:
                    res = copy.copy(atom.residue)
                    res.molecule = None
                    res.chain = atom.chain
                    res.children = ChildList(res)

                    res.chain.add(res)
                    replaced[atom.residue] = res
                else:
                    res = replaced[atom.residue]

                atom.residue = res
                assert atom.residue.chain is atom.chain
                res.add(atom)

            for oldnbr, bondorder in oldatom.bond_graph.items():
                if oldnbr not in old_to_new: continue
                newnbr = old_to_new[oldnbr]
                temp_bond_graph[atom][newnbr] = bondorder

        # Finally, do a deepcopy, which bring all the information along without linking it
        newatoms, new_bond_graph = copy.deepcopy((tempatoms, temp_bond_graph))

        for atom, original in zip(newatoms, oldatoms):
            atom.bond_graph = new_bond_graph[atom]
            atom.position = original.position
            atom.momentum = original.momentum

        return AtomList(newatoms)

    ###########################################
    # Routines to modify the geometry
    def rotate(self, angle, axis, center=None):
        """Rotate this object in 3D space

        Args:
            angle (u.Scalar[angle]): angle to rotate by
            axis (u.Vector[length]): axis to rotate about (len=3)
            center (u.Vector[length]): center of rotation (len=3) (default: origin)
        """
        center = utils.if_not_none(center, self.com)

        if hasattr(angle, 'units'): angle = angle.value_in(u.radians)
        rotmat = external.transformations.rotation_matrix(angle, axis, point=center)
        self.transform(rotmat)

    def translate(self, vector):
        """Translate this object in 3D space

        Args:
            vector (u.Vector[length]): translation vector, len=3
        """
        for atom in self.atoms:
            atom.position += vector

    def transform(self, matrix):
        """ Transform this object's coordinates using the provided 4x4 matrix

        Args:
            matrix (numpy.ndarray): transformation matrix, shape=(4,4)
        """
        # TODO: deal with units ... hard because the matrix has diff units for diff columns
        assert matrix.shape == (4, 4)
        self.positions = mathutils.apply_4x4_transform(matrix, self.positions)

    def atoms_within(self, radius, other=None, include_self=False):
        """ Return all atoms in an object within a given radius of this object

        Args:
            radius (u.Scalar[length]): radius to search for atoms
            other (AtomContainer): object containing the atoms to search (default:self.molecule)
            include_self (bool): if True, include the atoms from this object (since, by definition,
               their distance from this object is 0)

        Returns:
            AtomList: list of the atoms within ``radius`` of this object
        """
        if other is None:
            other = self.atoms[0].molecule
        if not include_self:
            filter_atoms = set(self.atoms)
        else:
            filter_atoms = set()

        distances = self.calc_distance_array(other=other)
        mindists = distances.min(axis=0)

        close_atoms = AtomList(atom for dist,atom in zip(mindists, other.atoms)
                               if dist <= radius and atom not in filter_atoms)
        return close_atoms

    def residues_within(self, radius, other=None, include_self=False):
        """ Return all atoms in an object within a given radius of this object

        Args:
            radius (u.Scalar[length]): radius to search for atoms
            other (AtomContainer): object containing the atoms to search (default:self.molecule)
            include_self (bool): if True, include the atoms from this object (since, by definition,
               their distance from this object is 0)

        Returns:
            AtomList: list of the atoms within ``radius`` of this object
        """
        atoms = self.atoms_within(radius, other=other, include_self=include_self)
        residues = collections.OrderedDict((atom.residue, None) for atom in atoms)
        return list(residues.keys())


class _AtomArray(object):
    def __init__(self, attrname):
        self.attrname = attrname

    def __get__(self, instance, owner):
        return u.array([getattr(atom, self.attrname) for atom in instance.atoms])

    def __set__(self, instance, value):
        assert len(value) == instance.num_atoms
        for atom, atomval in zip(instance.atoms, value):
            setattr(atom, self.attrname, atomval)


class AtomContainer(AtomGroup):
    """
    Mixin functions for NON-MOLECULE objects that have a list of atoms at``self.atoms``
    """
    positions = _AtomArray('position')
    masses = _AtomArray('mass')
    momenta = _AtomArray('momentum')
    velocities = _AtomArray('velocity')

    def __add__(self, other):
        l = mdt.AtomList(self.atoms)
        l.extend(other.atoms)
        return l

    @property
    def bond_graph(self):
        """ Dict[moldesign.Atom: List[moldesign.Atom]]: bond graph for all atoms in this object
        """
        return {atom: atom.bond_graph for atom in self.atoms}

    @property
    def bonds(self):
        """ Iterable[moldesign.Bond]: iterator over bonds from this object's atoms
        """
        bg = self.bond_graph
        for atom, nbrs in bg.items():
            for nbr, order in nbrs.items():
                if atom.index < nbr.index or nbr not in bg:
                    yield mdt.Bond(atom,nbr, order)

    def get_bond(self, a1, a2):
        return mdt.Bond(a1, a2, order=self.bond_graph[a1][a2])

    @property
    def internal_bonds(self):
        """ Iterable[moldesign.Bond]: iterator over bonds that connect two atoms in this object
        """
        bg = self.bond_graph
        for atom, nbrs in bg.items():
            for nbr, order in nbrs.items():
                if atom.index < nbr.index and nbr in bg:
                    yield mdt.Bond(atom, nbr, order)

    @property
    def external_bonds(self):
        """
        Iterable[moldesign.Bond]: iterator over bonds that bond these atoms to other atoms
        """
        bg = self.bond_graph
        for atom, nbrs in bg.items():
            for nbr, order in nbrs.items():
                if nbr not in bg:
                    yield mdt.Bond(atom, nbr, order)

    @property
    def bonded_atoms(self):
        """ List[moldesign.Atom]: list of external atoms this object is bonded to
        """
        bg = self.bond_graph
        atoms = []
        for atom, nbrs in bg.items():
            for nbr, order in nbrs.items():
                if nbr not in bg:
                    atoms.append(nbr)
        return atoms

    def bonds_to(self, other):
        """ Returns list of bonds between this object and another one

        Args:
            other (AtomContainer): other object

        Returns:
            List[moldesign.Bond]: bonds between this object and another
        """
        bonds = []
        otheratoms = set(other.atoms)
        for bond in self.internal_bonds:
            if bond.a1 in otheratoms or bond.a2 in otheratoms:
                bonds.append(bond)
        return bonds


@toplevel
class AtomList(AtomContainer, list):  # order is important, list will override methods otherwise
    """ A list of atoms with various helpful methods for creating and manipulating atom selections

    Args:
        atomlist (List[AtomContainer]): list of objects that are either atoms or contain a list of
           atoms at ``atomlist.atoms``
    """
    def __init__(self, atomlist=()):
        atoms = []
        for obj in atomlist:
            if hasattr(obj, 'atoms'):
                atoms.extend(obj.atoms)
            else:
                atoms.append(obj)
        super().__init__(atoms)

    def __getitem__(self, item):
        result = super().__getitem__(item)
        if isinstance(item, slice):
            return type(self)(result)
        else:
            return result

    def __getslice__(self, i, j):
        result = super().__getslice__(i, j)
        return type(self)(result)

    def __str__(self):
        return '[Atoms: %s]' % ', '.join(atom._shortstr() for atom in self)

    def __repr__(self):
        try:
            return '<AtomList: [%s]>' % ', '.join(atom._shortstr() for atom in self)
        except (KeyError, AttributeError):
            return '<AtomList at %x (__repr__ failed)>' % id(self)

    copy = AtomContainer.copy_atoms

    def intersection(self, *otherlists):
        """ Return a list of atoms that appear in all lists (including this one).

        Args:
            *otheriters (Iterable): one or more lists of atoms

        Returns:
            moldesign.AtomList: intersection of this lists with all passed lists. Preserves order
              in this list
        """
        s = set(self).intersection(*otherlists)
        return type(self)(o for o in self if o in s)

    def union(self, *otherlists):
        """ Return a list of atoms that appear in any lists (including this one).

        Args:
            *otherlists (Iterable): one or more lists of atoms

        Returns:
            moldesign.AtomList: union of this list of atoms with all passed lists of atoms.
               Equivalent to concatenating all lists then removing duplicates
        """
        found = set()
        newlist = type(self)()
        for item in itertools.chain(self, *otherlists):
            if item not in found:
                found.add(item)
                newlist.append(item)
        return newlist

    def unique(self):
        """ Return only unique atoms from this list

        Returns:
            moldesign.AtomList: copy of this list without any duplicates. Preserves order.
        """
        return self.union()

    def __sub__(self, other):
        otherset = set(other)
        return type(self)(atom for atom in self if atom not in otherset)

    # alias for self so that this works with AtomContainer methods
    @property
    def atoms(self):
        return self
