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
from scipy.spatial import distance as spd

import moldesign as mdt
import moldesign.utils
import moldesign.utils.callsigs
import moldesign.utils.classes
from moldesign import units as u, utils
from moldesign import data, utils, viewer, ui
from moldesign.external import transformations

__all__ = 'Atom AtomList Bond'.split()

class AtomContainer(object):
    """ Mixin functions for objects that have a ``self.atoms`` attribute with a list of atoms

    Attributes:
        atoms (List[Atom]): a list of atoms
    """
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
        other = utils.if_not_none(other, self)
        try:
            other_positions = other.atoms.position.value_in(u.ang)
        except AttributeError:
            other_positions = np.array([other.position.value_in(u.ang)])

        distances = spd.cdist(self.atoms.position.value_in(u.ang), other_positions)
        return distances * u.ang

    def calc_displacements(self):
        """ Calculate an array of displacements between all atoms in this object

        Returns:
            u.Array[length]: array of pairwise displacements between atoms

        Example:
            >>> displacements = self.calc_displacements(other)
            >>> displacements[i, j] == (self.atoms[i].position - self.atoms[j].position)
        """
        # TODO: allow other, similar to calc_distance array
        return utils.pairwise_displacements(self.atoms.position)

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
            other = self.atoms[0].parent
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
    natoms = numatoms = num_atoms

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

    def draw(self, width=500, height=500, show_2dhydrogens=False):
        """ Visualize this molecule (Jupyter only).

        Creates a 3D viewer, and, for small molecules, a 2D viewer).

        Args:
            width (int): width of the viewer in pixels
            height (int): height of the viewer in pixels

        Returns:
            moldesign.ui.SelectionGroup
        """
        import ipywidgets as ipy
        viz2d = self.draw2d(width=width, height=height,
                            display=False,
                            show_hydrogens=show_2dhydrogens)
        viz3d = self.draw3d(width=width, height=height,
                            display=False)
        views = ipy.HBox([viz2d, viz3d])
        return ui.SelectionGroup([views, ui.AtomInspector()])

    def draw3d(self, highlight_atoms=None, **kwargs):
        """ Draw this object in 3D. Jupyter only.

        Args:
            highlight_atoms (List[Atom]): atoms to highlight when the structure is drawn

        Returns:
            mdt.GeometryViewer: 3D viewer object
        """
        self.viz3d = viewer.GeometryViewer(self, **kwargs)
        if highlight_atoms is not None:
            self.viz3d.highlight_atoms(highlight_atoms)
        return self.viz3d

    def draw2d(self, highlight_atoms=None, show_hydrogens=False, **kwargs):
        """
        Draw this object in 2D. Jupyter only.

        Args:
            highlight_atoms (List[Atom]): atoms to highlight when the structure is drawn

        Returns:
            mdt.ChemicalGraphViewer: 2D viewer object
        """
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
            atom.parent = None
            atom.bond_graph = {}

            if atom.chain is not None:
                if atom.chain not in replaced:
                    chain = copy.copy(atom.chain)
                    chain.parent = None
                    chain.children = {}
                    replaced[atom.chain] = chain
                else:
                    chain = replaced[atom.chain]

                atom.chain = chain

            if atom.residue is not None:
                if atom.residue not in replaced:
                    res = copy.copy(atom.residue)
                    res.parent = None
                    res.chain = atom.chain
                    res.children = {}

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
        rotmat = transformations.rotation_matrix(angle, axis, point=center)
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
        positions = u.to_units_array(self.atoms.position)
        newpos = mdt.geometry._apply_4trans(matrix, positions)
        for atom, r in zip(self.atoms, newpos):
            atom.position = r


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
        if a1.index > a2.index:
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
        try: ff = self.a1.parent.energy_model.get_forcefield()
        except (NotImplementedError, AttributeError): return None
        return ff.bond_term[self]

    def __repr__(self):
        try:
            return '<Bond: %s>'%str(self)
        except:
            print '<Bond @ %s (error in __repr__)>' % id(self)

    def __str__(self):
        return self.name


class ProtectedArray(object):
    """
    Descriptor for arrays that shouldn't be reassigned.
    Makes sure array attributes (specifically position and momentum) are modified in place

    Args:
        name (str): name of the instance attribute
    """
    def __init__(self, name):
        self.name = name

    def __get__(self, instance, cls=None):
        return getattr(instance, self.name)

    def __set__(self, instance, value):
        self.__get__(instance)[:] = value


class AtomArray(ProtectedArray):
    """
    Descriptor for atom coordinates are stored in the parent molecule.

    Makes sure that the arrays and their references are maintained during both
    reassignment and copying/pickling

    Args:
        atomname (str): name of the attribute in the atom instance
        parentname (str): name of the corresponding attribute in the molecule instance
    """
    def __init__(self, atomname, parentname):
        self.name = atomname
        self.parentname = parentname

    def __get__(self, instance, cls=None):
        if instance.parent is None:
            return getattr(instance, self.name)
        else:
            return getattr(instance.parent, self.parentname)[instance.parent_slice]


class AtomCoordinate(object):
    """ Descriptor for use with the ``Atom`` class.

    Gives access to 3D coordinates as ``atom.x,atom.y,atom.z`` instead of
    ``atom.position[0],atom.position[1],atom.position[2]``

    Args:
        quantity (str): name of the attribute that this accesses
        index (int): component of the attribute that this accesses
    """
    def __init__(self, attrname, index):
        self.attrname = attrname
        self.index = index

    def __get__(self, instance, cls):
        array = getattr(instance, self.attrname)
        return array[self.index]

    def __set__(self, instance, value):
        array = getattr(instance, self.attrname)
        array[self.index] = value


class AtomDrawingMixin(object):
    """ Functions for creating atomic visualizations.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Atom` class. Routines
        are separated are here for code organization only - they could be included in the main
        Atom class without changing any functionality
    """

    #@utils.args_from(mdt.molecule.Molecule.draw2d, allexcept=['highlight_atoms'])  # import order
    def draw2d(self, **kwargs):
        """ Draw a 2D viewer with this atom highlighted (Jupyter only).
        In biomolecules, only draws the atom's residue.

        Args:
            width (int): width of viewer in pixels
            height (int): height of viewer in pixels

        Returns:
            mdt.ChemicalGraphViewer: viewer object
        """
        if self.parent:
            if self.parent.is_small_molecule:
                return self.parent.draw2d(highlight_atoms=[self], **kwargs)
            elif self.parent.is_biomolecule:
                return self.residue.draw2d(highlight_atoms=[self], **kwargs)
            else:
                raise ValueError('No drawing routine specified')
        else:
            raise ValueError('No drawing routine specified')

    #@utils.args_from(mdt.molecule.Molecule.draw2d, allexcept=['highlight_atoms'])  # import order
    def draw3d(self, **kwargs):
        """ Draw a 3D viewer with this atom highlighted (Jupyter only).

        Args:
            width (int): width of viewer in pixels
            height (int): height of viewer in pixels

        Returns:
            mdt.GeometryViewer: viewer object
        """
        return self.parent.draw3d(highlight_atoms=[self], **kwargs)

    def draw(self, width=300, height=300):
        """ Draw a 2D and 3D viewer with this atom highlighted (notebook only)

        Args:
            width (int): width of viewer in pixels
            height (int): height of viewer in pixels

        Returns:
            ipy.HBox: viewer object
        """
        import ipywidgets as ipy
        viz2d = self.draw2d(width=width, height=height, display=False)
        viz3d = self.draw3d(width=width, height=height, display=False)
        return ipy.HBox([viz2d, viz3d])


class AtomGeometryMixin(object):
    """ Functions measuring distances between atoms and other things.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Atom` class. Routines
        are separated are here for code organization only - they could be included in the main
        Atom class without changing any functionality
    """
    @moldesign.utils.callsigs.args_from(AtomContainer.distance)
    def distance(self, *args, **kwargs):
        return self._container.distance(*args, **kwargs)

    @moldesign.utils.callsigs.args_from(AtomContainer.atoms_within)
    def atoms_within(self, *args, **kwargs):
        return self._container.atoms_within(*args, **kwargs)

    @moldesign.utils.callsigs.args_from(AtomContainer.calc_distance_array)
    def calc_distances(self, *args, **kwargs):
        array = self._container.calc_distance_array(*args, **kwargs)
        return array[0]

    @property
    def _container(self):
        """ AtomContainer: a container with just this atom in it.

        This is a convenience method for accessing to all of the :class:`AtomContainer`'s
        useful methods for dealing with geometry
        """
        return AtomList([self])


class AtomPropertyMixin(object):
    """ Functions accessing computed atomic properties.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Atom` class. Routines
        are separated are here for code organization only - they could be included in the main
        Atom class without changing any functionality
    """
    @property
    def ff(self):
        """ moldesign.utils.DotDict: This atom's force field parameters, if available (``None``
        otherwise)
        """
        try:
            ff = self.parent.energy_model.mdtforcefield
        except AttributeError:
            return None
        if ff is None: return None

        return moldesign.utils.classes.DotDict(partialcharge=ff.partial_charges[self],
                                               lj=ff.lennard_jones[self])

    @property
    def basis_functions(self):
        """ List[mdt.orbitals.AtomicBasisFunction]: This atom's basis functions, if available (
        ``None`` otherwise)
        """
        if self.parent is None:
            return None
        try:
            wfn = self.parent.electronic_state
        except mdt.molecule.NotCalculatedError:
            return None

        return wfn.aobasis.on_atom.get(self, [])

    @property
    def properties(self):
        """ moldesign.utils.DotDict: Returns any calculated properties for this atom
        """
        props = moldesign.utils.classes.DotDict()
        for name, p in self.parent.properties.iteritems():
            if hasattr(p, 'type') and p.type == 'atomic':
                props[name] = p[self]
        return props


class AtomReprMixin(object):
    """ Functions for printing out various strings related to the atom.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Atom` class. Routines
        are separated are here for code organization only - they could be included in the main
        Atom class without changing any functionality
    """
    def __str__(self):
        desc = '%s %s (elem %s)' % (self.__class__.__name__, self.name, self.elem)
        molstring = ''
        if self.parent:
            molstring = ', index %d' % self.index
            if self.parent.is_biomolecule:
                molstring += ' (res %s chain %s)' % (self.residue.name, self.chain.name)
        return '%s%s' % (desc, molstring)

    def __repr__(self):
        # TODO: rename parent to "molecule"
        try:
            if self.parent:
                return '<%s in molecule %s>' % (self, self.parent)
            else:
                return '<%s>' % self
        except:
            return '<%s at %s (exception in __repr__)>' % (self.__class__.__name__, id(self))

    def markdown_summary(self):
        """Return a markdown-formatted string describing this atom

        Returns:
            str: markdown-formatted string
            """
        if self.parent is None:
            lines = ["<h3>Atom %s</h3>" % self.name]
        else:
            lines = ["<h3>Atom %s (index %d)</h3>" % (self.name, self.index)]

        lines.append('**Atomic number**: %d' % self.atnum)
        lines.append("**Mass**: %s" % self.mass)

        if self.parent is not None:
            if self.parent.is_biomolecule:
                lines.append("\n\n**Residue**: %s (index %d)" % (self.residue.name, self.residue.index))
                lines.append("**Chain**: %s" % self.chain.name)
            lines.append("**Molecule**: %s" % self.parent.name)
            for ibond, (nbr, order) in enumerate(self.bond_graph.iteritems()):
                lines.append('**Bond %d** (order = %d): %s (index %s) in %s' % (
                    ibond + 1, order, nbr.name, nbr.index, nbr.residue.name))

        if self.basis_functions:
            lines.append('**Basis functions:**<br>' + '<br>'.join(map(str,self.basis_functions)))

        if self.ff:
            lines.append('**Forcefield partial charge**: %s' % self.ff.partialcharge)
            # TODO: deal with other LJ types, e.g., AB?
            lines.append(u'**Forcefield LJ params**: '
                         u'\u03C3=%s, \u03B5=%s' % (
                             self.ff.lj.sigma.defunits(),
                             self.ff.lj.epsilon.defunits()))

        # position and momentum
        if self.ndim == 3:
            table = utils.MarkdownTable('', 'x', 'y', 'z')
        else:
            table = utils.MarkdownTable(*([''] + map(str, range(self.ndim))))
        table.add_line(['**position /** {}'.format(u.default.length)] +
                       ['%12.3f' % x.defunits_value() for x in self.position])
        table.add_line(['**momentum /** {}'.format(u.default.momentum)] +
                       ['%12.3e' % m.defunits_value() for m in self.momentum])
        try:
            self.force
        except:
            pass
        else:
            table.add_line(['**force /** {.units}'.format(self.force.defunits())] +
                           ['%12.3e' % m.defunits_value() for m in self.force])

        lines.append('\n\n' + table.markdown() + '\n\n')
        # All other assigned properties

        return '<br>'.join(lines)

    def _repr_markdown_(self):
        return self.markdown_summary()


class Atom(AtomDrawingMixin, AtomGeometryMixin, AtomPropertyMixin, AtomReprMixin):
    """
    Holds atomic info.
    Once assigned to a Molecule, position and momentum are automatically linked to the molecule's
    coordinates.

    Args:
        name (str): The atom's name (if not passed, set to the element name + the atom's index)
        atnum (int): Atomic number (if not passed, determined from element if possible)
        mass (u.Scalar[mass]): The atomic mass (if not passed, set to the most abundant isotopic
            mass)
        chain (moldesign.Chain): biomolecular chain that this atom belongs to
        residue (moldesign.Residue): biomolecular residue that this atom belongs to
        pdbname (str): name from PDB entry, if applicable
        pdbindex (int): atom serial number in the PDB entry, if applicable
        element (str): Elemental symbol (if not passed, determined from atnum if possible)

    Attributes:
        name (str): A descriptive name for this atom
        element (str): IUPAC elemental symbol ('C', 'H', 'Cl', etc.)
        index (int): the atom's current index in the molecule
            (``self is self.parent.atoms[ self.index]``)

        position (u.Vector[length, 3]): atomic positions
        momentum (u.Vector[momentum, 3]): atomic momenta
        atnum (int): atomic number (synonyms: atomic_num)
        mass (u.Scalar[mass]): the atom's mass

        residue (moldesign.Residue): biomolecular residue that this atom belongs to
        chain (moldesign.Chain): biomolecular chain that this atom belongs to
        parent (moldesign.Molecule): molecule that this atom belongs to
        index (int): index in the parent molecule: ``atom is atom.parent.atoms[index]``

        x,y,z (u.Scalar[length]): x, y, and z components of ``atom.position``
        vx, vy, vz (u.Scalar[length/time]): x, y, of ``atom.velocity``
        px, py, pz (u.Scalar[momentum]): x, y, and z of ``atom.momentum``
        fx, fy, fz (u.Scalar[force]): x, y, and z ``atom.force``
            Raises:
                moldesign.molecule.NotCalculatedError: if forces are not calculated
    """
    x, y, z = (AtomCoordinate('position', i) for i in xrange(3))
    vx, vy, vz = (AtomCoordinate('velocity', i) for i in xrange(3))
    px, py, pz = (AtomCoordinate('momentum', i) for i in xrange(3))
    fx, fy, fz = (AtomCoordinate('force', i) for i in xrange(3))
    position = AtomArray('_position', 'positions')
    momentum = AtomArray('_momentum', 'momenta')

    #################################################################
    # Methods for BUILDING the atom and indexing it in a molecule
    def __init__(self, name=None, atnum=None, mass=None, chain=None, residue=None,
                 pdbname=None, pdbindex=None, element=None):
        if element: self.atnum = data.ATOMIC_NUMBERS[element]
        else: self.atnum = atnum

        self.name = utils.if_not_none(name, self.elem)
        self.pdbname = utils.if_not_none(pdbname, self.name)
        self.pdbindex = pdbindex

        if mass is None: self.mass = data.ATOMIC_MASSES[self.atnum]
        else: self.mass = mass

        self.residue = residue
        self.chain = chain
        self.parent = None
        self.index = None
        self._parent_slice = None
        self._position = np.zeros(self.ndim) * u.default.length
        self._momentum = np.zeros(self.ndim) * (u.default.length*
                                                u.default.mass/u.default.time)
        self._bond_graph = {}

    @moldesign.utils.callsigs.args_from(AtomContainer.copy)
    def copy(self, *args, **kwargs):
        return self._container.copy(*args, **kwargs)[0]

    def __getstate__(self):
        state = self.__dict__.copy()
        if self.parent is not None:  # these don't belong to the atom anymore
            state['_bond_graph'] = None
            state['_position'] = self.position
            state['_momentum'] = self.momentum
        return state

    def _set_parent(self, parent):
        """ Permanently bind this atom to a parent molecule

        Args:
            parent (moldesign.Molecule): the molecule that this atom will become a part of
        """
        if self.parent and (parent is not self.parent):
            raise ValueError('Atom is already part of a molecule')
        self.parent = parent

    def _index_into_molecule(self, array_name, moleculearray, molslice):
        """ Change the atom's position/momentum etc. arrays to be pointers into
        the appropriate slices of the parent molecule's master array of positions/momenta

        Args:
            array_name: the private name of the array (assumes private name is '_'+array_name)
            moleculearray: the molecule's master array
            molslice: the python slice object

        Note:
            This will be called by the molecule's init method
        """
        oldarray = getattr(self, array_name)
        moleculearray[molslice] = oldarray
        setattr(self, '_' + array_name, None)

    def bond_to(self, other, order):
        """ Create a bond to another atom

        Args:
            other (Atom): atom to bond to
            order (int): bond order

        Returns:
            Bond: bond object
        """
        if self.parent is other.parent:
            self.bond_graph[other] = other.bond_graph[self] = order
            if self.parent is not None: self.parent.num_bonds += 1
        else:  # allow unassigned atoms to be bonded to anything for building purposes
            assert self.parent is None
            self.bond_graph[other] = order
        return Bond(self, other, order)

    @property
    def bond_graph(self):
        """ Mapping[Atom, int]: a dictionary containing all other atoms this atom is bonded to,
        of the form
           ``{nbr1: bond_order_1, nbr2: bond_order_2, ...}``
        """
        if self.parent is None:
            return self._bond_graph
        else:
            self._bond_graph = None
            try:
                return self.parent.bond_graph[self]
            except KeyError:
                self.parent.bond_graph[self] = {}
                return self.parent.bond_graph[self]

    @bond_graph.setter
    def bond_graph(self, value):
        if self.parent is None:
            self._bond_graph = value
        else:
            self._bond_graph = None
            self.parent.bond_graph[self] = value

    @property
    def bonds(self):
        """ List[Bond]: list of all bonds this atom is involved in
        """
        return [Bond(self, nbr, order) for nbr,order in self.bond_graph.iteritems()]

    @property
    def force(self):
        """ u.Vector[force, 3]: force on this atom
        """
        f = self.parent.forces
        return f[self._parent_slice]

    @property
    def velocity(self):
        """ u.Vector[length/time, 3]: velocity of this atom
        """
        return (self.momentum / self.mass).defunits()

    @property
    def num_bonds(self):
        """ int: the number of other atoms this atom is bonded to
        """
        return len(self.bond_graph)

    @property
    def symbol(self):
        """ str: elemental symbol
        """
        return data.ELEMENTS.get(self.atnum,'?')
    elem = element = symbol


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
            return u.to_units_array(atom_attrs)
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


# Hack for documentation only - sure these get picked up in the docs
if mdt._building_docs:
    __all__.extend(('AtomDrawingMixin AtomGeometryMixin AtomPropertyMixin AtomReprMixin '
                    'AtomContainer').split())
