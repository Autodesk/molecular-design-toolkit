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

import copy

import numpy as np

import moldesign as mdt
from moldesign import units as u
from moldesign import utils, external, mathutils
from . import toplevel


class AtomContainer(object):
    """ Mixin functions for objects that have a ``self.atoms`` attribute with a list of atoms

    Attributes:
        atoms (List[Atom]): a list of atoms
    """

    @property
    def positions(self):
        """ u.Array[length]: (Nx3) array of atomic positions
        """
        return self.atoms.position

    @positions.setter
    def positions(self, val):
        assert len(val) == self.num_atoms
        for atom, p in zip(self.atoms, val):
            atom.position = p

    def __init__(self, *args, **kwargs):
        """ This should never be called directly - it will be called by the `super` methods
        of its subclasses """
        super(AtomContainer, self).__init__(*args, **kwargs)
        self._atom_attrs = None
        self.viz2d = None
        self.viz3d = None

    @property
    def heavy_atoms(self):
        """ AtomList: a list of all heavy atoms (i.e., non-hydrogen) in this object """
        return AtomList([a for a in self.atoms if a.atnum != 1])

    @property
    def mass(self):
        """ u.Scalar[mass]: total mass of this object
        """
        return sum(a.mass for a in self.atoms)

    def get_atoms(self, **queries):
        """Allows keyword-based atom queries.

        Args:
            **queries (dict): parameters to match

        Returns:
            AtomList: the atoms matching this query
        """
        if not queries: return self

        def atom_callback(atom, attrs):
            obj = atom
            for attr in attrs:
                obj = getattr(obj, attr)
            return obj

        result = AtomList()
        for atom in self.atoms:
            for q, matches in queries.iteritems():
                attr = atom_callback(atom, q.split('.'))
                if type(matches) in (list, tuple):
                    if attr not in matches: break
                else:
                    if attr != matches: break
            else:  # if here, this atom matched all queries
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

    def atoms_within(self, radius, other=None, include_self=False):
        """ Return all atoms in an object within a given radius of this object

        Args:
            radius (u.Scalar[length]): radius to search for atoms
            other (AtomContainer): object containing the atoms to search (default:self.parent)
            include_self (bool): if True, include the atoms from this object (since, by definition,
               their distance from this object is 0)

        Returns:
            AtomList: list of the atoms within ``radius`` of this object
        """
        if other is None:
            other = self.atoms[0].molecule
        if not include_self:
            myatoms = set(self.atoms)

        close_atoms = AtomList()
        for atom in other.atoms:
            if self.distance(atom) <= radius and (include_self or (atom not in myatoms)):
                close_atoms.append(atom)
        return close_atoms

    @property
    def num_atoms(self):
        """ int: number of atoms in this object """
        return len(self.atoms)
    natoms = num_atoms

    @property
    def center_of_mass(self):
        """ units.Vector[length]: The (x,y,z) coordinates of this object's center of mass """
        total_mass = 0.0 * u.default.mass
        com = np.zeros(3) * u.default.length * u.default.mass
        for atom in self.atoms:
            total_mass += atom.mass
            com += atom.position * atom.mass
        com = com / total_mass
        return com
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
        return mdt.geom.angle(*map(self._getatom, (a1, a2, a3)))

    def dihedral(self, a1, a2, a3=None, a4=None):
        """ Calculate the dihedral angle between atoms a1, a2, a3, a4.

        Atoms can be passed as the atoms themselves or as the atom names

        Args:
            a1, a2, a3, a4 (str OR int OR moldesign.Atom): atoms defining the dihedral

        Returns:
            units.Scalar[angle]
        """
        return mdt.geom.dihedral(*map(self._getatom, (a1, a2, a3, a4)))

    def draw(self, width=500, height=500, show_2dhydrogens=None, display=False):
        """ Visualize this molecule (Jupyter only).

        Creates a 3D viewer, and, for small molecules, a 2D viewer).

        Args:
            width (int): width of the viewer in pixels
            height (int): height of the viewer in pixels
            show_2dhydrogens (bool): whether to show the hydrogens in 2d (default: True if there
                   are 10 or less heavy atoms, false otherwise)
            display (bool): immediately display this viewer

        Returns:
            moldesign.ui.SelectionGroup
        """
        import ipywidgets as ipy
        import IPython.display

        if self.num_atoms < 40:

            viz2d = self.draw2d(width=width, height=height,
                                display=False,
                                show_hydrogens=show_2dhydrogens)
            viz3d = self.draw3d(width=width, height=height,
                                display=False)
            views = ipy.HBox([viz2d, viz3d])
        else:
            views = self.draw3d(display=False)

        displayobj = mdt.uibase.SelectionGroup([views, mdt.uibase.components.AtomInspector()])

        if display:
            IPython.display.display(displayobj)
        return displayobj

    def draw3d(self, highlight_atoms=None, **kwargs):
        """ Draw this object in 3D. Jupyter only.

        Args:
            highlight_atoms (List[Atom]): atoms to highlight when the structure is drawn

        Returns:
            mdt.GeometryViewer: 3D viewer object
        """
        from moldesign import viewer
        self.viz3d = viewer.GeometryViewer(self, **kwargs)
        if highlight_atoms is not None:
            self.viz3d.highlight_atoms(highlight_atoms)
        return self.viz3d

    def draw2d(self, highlight_atoms=None, show_hydrogens=None, **kwargs):
        """
        Draw this object in 2D. Jupyter only.

        Args:
            highlight_atoms (List[Atom]): atoms to highlight when the structure is drawn
            show_hydrogens (bool): whether to draw the hydrogens or not (default: True if there
                   are 10 or less heavy atoms, false otherwise)

        Returns:
            mdt.ChemicalGraphViewer: 2D viewer object
        """
        from moldesign import viewer
        if show_hydrogens is None:
            show_hydrogens = len(self.heavy_atoms) <= 10
        if not show_hydrogens:
            alist = [atom for atom in self.atoms if atom.atnum > 1]
        else:
            alist = self
        self.viz2d = viewer.DistanceGraphViewer(alist, **kwargs)
        if highlight_atoms: self.viz2d.highlight_atoms(highlight_atoms)
        return self.viz2d

    def copy(self):
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

            for oldnbr, bondorder in oldatom.bond_graph.iteritems():
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
        positions = u.array(self.atoms.position)
        newpos = mathutils.apply_4x4_transform(matrix, positions)
        for atom, r in zip(self.atoms, newpos):
            atom.position = r


@toplevel
class AtomList(list, AtomContainer):
    """A list of atoms that allows attribute "slicing" - accessing an attribute of the
    list will return a list of atom attributes.

    Example:
        >>> atomlist.mass == [atom.mass for atom in atomlist.atoms]
        >>> getattr(atomlist, attr) == [getattr(atom, attr) for atom in atomlist.atoms]
    """

    def __init__(self, *args, **kwargs):
        super(AtomList, self).__init__(*args, **kwargs)

    def __dir__(self):
        """
        Overrides ``dir`` to allow tab completion for both this object's attributes
        and the underlying atomic properties
        """
        if self._atom_attrs is None:
            self._atom_attrs = set([k for k,v in self.atoms[0].__dict__.iteritems()
                                    if not callable(v)])
        return list(self._atom_attrs.union(dir(self.__class__)).union(self.__dict__))

    def __getattr__(self, item):
        """
        Returns an array of atomic properties, e.g., an array of all atomic masses
        is returned by ``AtomContainer.mass = AtomContainer.__getattr__('mass')``
        """
        if item.startswith('__'):
            raise AttributeError
        atom_attrs = [getattr(atom, item) for atom in self.atoms]
        try:
            return u.array(atom_attrs)
        except (TypeError, StopIteration):
            return atom_attrs

    def __setattr__(self, key, vals):
        """
        Set an array of atomic properties, e.g., set all atomic masses to 1 with
        ``atoms.mass = [1.0 for a in atoms]``
        """
        assert len(vals) == len(self)
        for atom, v in zip(self, vals): setattr(atom, key, v)

    @property
    def atoms(self):
        """This is a synonym for self so that AtomContainer methods will work here too"""
        return self

    def __getslice__(self, *args, **kwargs):
        return AtomList(super(AtomList, self).__getslice__(*args, **kwargs))