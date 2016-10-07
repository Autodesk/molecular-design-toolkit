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

import numpy as np

import moldesign as mdt
from moldesign import data, utils
from moldesign import units as u
from . import toplevel, AtomContainer, AtomList, AtomArray, AtomCoordinate, Bond


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
        if self.molecule:
            if self.molecule.is_small_molecule:
                return self.molecule.draw2d(highlight_atoms=[self], **kwargs)
            elif self.molecule.is_biomolecule:
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
        return self.molecule.draw3d(highlight_atoms=[self], **kwargs)

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
    @utils.args_from(AtomContainer.distance)
    def distance(self, *args, **kwargs):
        return self._container.distance(*args, **kwargs)

    @utils.args_from(AtomContainer.atoms_within)
    def atoms_within(self, *args, **kwargs):
        return self._container.atoms_within(*args, **kwargs)

    @utils.args_from(AtomContainer.calc_distance_array)
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
            ff = self.molecule.energy_model.mdtforcefield
        except AttributeError:
            return None
        if ff is None: return None

        return utils.DotDict(partialcharge=ff.partial_charges[self],
                             lj=ff.lennard_jones[self])

    @property
    def basis_functions(self):
        """ List[mdt.orbitals.AtomicBasisFunction]: This atom's basis functions, if available
        (``None`` otherwise)
        """
        if self.molecule is None:
            return None
        try:
            wfn = self.molecule.wfn
        except mdt.exceptions.NotCalculatedError:
            return None

        return wfn.aobasis.on_atom.get(self, [])

    @property
    def properties(self):
        """ moldesign.utils.DotDict: Returns any calculated properties for this atom
        """
        props = utils.DotDict()
        for name, p in self.molecule.properties.iteritems():
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
        if self.molecule:
            molstring = ', index %d' % self.index
            if self.molecule.is_biomolecule:
                molstring += ' (res %s chain %s)' % (self.residue.name, self.chain.name)
        return '%s%s' % (desc, molstring)

    def __repr__(self):
        # TODO: rename parent to "molecule"
        try:
            if self.molecule:
                return '<%s in molecule %s>' % (self, self.molecule)
            else:
                return '<%s>' % self
        except:
            return '<%s at %s (exception in __repr__)>' % (self.__class__.__name__, id(self))

    def markdown_summary(self):
        """Return a markdown-formatted string describing this atom

        Returns:
            str: markdown-formatted string
        """
        if self.molecule is None:
            lines = ["<h3>Atom %s</h3>" % self.name]
        else:
            lines = ["<h3>Atom %s (index %d)</h3>" % (self.name, self.index)]

        lines.append('**Atomic number**: %d' % self.atnum)
        lines.append("**Mass**: %s" % self.mass)
        lines.append('**Formal charge**: %s' % self.formal_charge)

        if self.molecule is not None:
            if self.molecule.is_biomolecule:
                lines.append("\n\n**Residue**: %s (index %d)" % (self.residue.name, self.residue.index))
                lines.append("**Chain**: %s" % self.chain.name)
            lines.append("**Molecule**: %s" % self.molecule.name)
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
        table = utils.MarkdownTable('', 'x', 'y', 'z')

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


@toplevel
class Atom(AtomDrawingMixin, AtomGeometryMixin, AtomPropertyMixin, AtomReprMixin):
    """ A data structure representing an atom.

    ``Atom`` objects store information about individual atoms within a larger molecular system,
    providing access to atom-specific geometric, biomolecular, topological and property
    information. Each :class:`Molecule<moldesign.Molecule>` is composed of a list of atoms.

    Atoms can be instantiated directly, but they will generally be created
    automatically as part of molecules.

    Args:
        name (str): The atom's name (if not passed, set to the element name + the atom's index)
        atnum (int): Atomic number (if not passed, determined from element if possible)
        mass (units.Scalar[mass]): The atomic mass (if not passed, set to the most abundant isotopic
            mass)
        chain (moldesign.Chain): biomolecular chain that this atom belongs to
        residue (moldesign.Residue): biomolecular residue that this atom belongs to
        pdbname (str): name from PDB entry, if applicable
        pdbindex (int): atom serial number in the PDB entry, if applicable
        element (str): Elemental symbol (if not passed, determined from atnum if possible)

    **Atom instance attributes:**

    Attributes:
        name (str): A descriptive name for this atom
        element (str): IUPAC elemental symbol ('C', 'H', 'Cl', etc.)
        index (int): the atom's current index in the molecule
            (``self is self.parent.atoms[ self.index]``)
        atnum (int): atomic number (synonyms: atomic_num)
        mass (u.Scalar[mass]): the atom's mass

        position (units.Vector[length]): atomic position vector. Once an atom is part of a molecule,
            this quantity will refer to ``self.molecule.positions[self.index]``.
        momentum (units.Vector[momentum]): atomic momentum vector. Once an atom is part of a
           molecule, this quantity will refer to ``self.molecule.momenta[self.index]``.

        x,y,z (u.Scalar[length]): x, y, and z components of ``atom.position``
        vx, vy, vz (u.Scalar[length/time]): x, y, of ``atom.velocity``
        px, py, pz (u.Scalar[momentum]): x, y, and z of ``atom.momentum``
        fx, fy, fz (u.Scalar[force]): x, y, and z ``atom.force``

        residue (moldesign.Residue): biomolecular residue that this atom belongs to
        chain (moldesign.Chain): biomolecular chain that this atom belongs to
        parent (moldesign.Molecule): molecule that this atom belongs to
        index (int): index in the parent molecule: ``atom is atom.parent.atoms[index]``

    **Atom methods and properties**

    See also methods offered by the mixin superclasses:

            - :class:`AtomDrawingMixin`
            - :class:`AtomGeometryMixin`
            - :class:`AtomPropertyMixin`
            - :class:`AtomReprMixin`
    """
    x, y, z = (AtomCoordinate('position', i) for i in xrange(3))
    vx, vy, vz = (AtomCoordinate('velocity', i) for i in xrange(3))
    px, py, pz = (AtomCoordinate('momentum', i) for i in xrange(3))
    fx, fy, fz = (AtomCoordinate('force', i) for i in xrange(3))
    position = AtomArray('_position', 'positions')
    momentum = AtomArray('_momentum', 'momenta')

    atomic_number = utils.Synonym('atnum')

    #################################################################
    # Methods for BUILDING the atom and indexing it in a molecule
    def __init__(self, name=None, atnum=None, mass=None, chain=None, residue=None,
                 formal_charge=None, pdbname=None, pdbindex=None, element=None):

        # Allow user to instantiate an atom as Atom(6) or Atom('C')
        if atnum is None and element is None:
            if isinstance(name, int):
                atnum = name
                name = None
            else: element = name

        if element: self.atnum = data.ATOMIC_NUMBERS[element]
        else: self.atnum = atnum

        self.name = utils.if_not_none(name, self.elem)
        self.pdbname = utils.if_not_none(pdbname, self.name)
        self.pdbindex = pdbindex

        if mass is None: self.mass = data.ATOMIC_MASSES[self.atnum]
        else: self.mass = mass

        self.formal_charge = utils.if_not_none(formal_charge, 0.0 * u.q_e)
        self.residue = residue
        self.chain = chain
        self.molecule = None
        self.index = None
        self._position = np.zeros(3) * u.default.length
        self._momentum = np.zeros(3) * (u.default.length*
                                       u.default.mass/u.default.time)
        self._bond_graph = {}

    def to_json(self, parent=None):
        """Designed to be called by the MdtJsonEncoder"""
        js = mdt.chemjson.jsonify(self,
                                  ('name index atnum position '
                                   'mass momentum').split())
        js['chain'] = self.chain.index
        js['residue'] = self.residue.index
        return js

    @utils.args_from(AtomContainer.copy)
    def copy(self, *args, **kwargs):
        """ Copy an atom (delegate to AtomContainer)
        """
        return self._container.copy(*args, **kwargs)[0]

    def __getstate__(self):
        """Helper for pickling"""
        state = self.__dict__.copy()
        if self.molecule is not None:  # then these don't belong to the atom anymore
            state['_bond_graph'] = None
            state['_position'] = self.position
            state['_momentum'] = self.momentum
        return state

    def _set_molecule(self, molecule):
        """ Permanently make this atom part of a molecule (private)

        Args:
            parent (moldesign.Molecule): the molecule that this atom will become a part of
        """
        if self.molecule and (molecule is not self.molecule):
            raise ValueError('Atom is already part of a molecule')
        self.molecule = molecule

    def _index_into_molecule(self, array_name, moleculearray, index):
        """ Link the atom's positions and momenta to the parent molecule (private)

        Args:
            array_name (str): the private name of the array (assumes private name is '_'+array_name)
            moleculearray (u.Array): the molecule's master array
            index: the atom's index in the molecule

        Note:
            This will be called by the molecule's init method
        """
        oldarray = getattr(self, array_name)
        moleculearray[index, :] = oldarray
        setattr(self, '_' + array_name, None)  # remove the internally stored version

    def bond_to(self, other, order):
        """ Create or modify a bond with another atom

        Args:
            other (Atom): atom to bond to
            order (int): bond order

        Returns:
            moldesign.molecules.bonds.Bond: bond object
        """
        if self.molecule is other.molecule:
            self.bond_graph[other] = other.bond_graph[self] = order
        else:  # allow unassigned atoms to be bonded to anything for building purposes
            self.bond_graph[other] = order
        return Bond(self, other, order)

    @property
    def bond_graph(self):
        """ Mapping[Atom, int]: dictionary of this atoms bonded neighbors, of the form
        ``{bonded_atom1, bond_order1, ...}``
        """
        if self.molecule is None:
            return self._bond_graph
        else:
            self._bond_graph = None
            try:
                return self.molecule.bond_graph[self]
            except KeyError:
                self.molecule.bond_graph[self] = {}
                return self.molecule.bond_graph[self]

    @bond_graph.setter
    def bond_graph(self, value):
        if self.molecule is None:
            self._bond_graph = value
        else:
            self._bond_graph = None
            self.molecule.bond_graph[self] = value

    @property
    def bonds(self):
        """ List[Bond]: list of all bonds this atom is involved in
        """
        return [Bond(self, nbr, order) for nbr, order in self.bond_graph.iteritems()]

    @property
    def heavy_bonds(self):
        """ List[Bond]: list of all heavy atom bonds (where BOTH atoms are not hydrogen)

        Note:
            this returns an empty list if called on a hydrogen atom
        """
        if self.atnum == 1:
            return []
        else:
            return [Bond(self, nbr, order)
                    for nbr, order in self.bond_graph.iteritems()
                    if nbr.atnum > 1]

    @property
    def force(self):
        """ (units.Vector[force]): atomic force vector. This quantity must be calculated - it is
        equivalent to ``self.molecule.forces[self.index]``

        Raises:
            moldesign.NotCalculatedError: if molecular forces have not been calculated
        """
        f = self.molecule.forces
        return f[self.index]

    @property
    def velocity(self):
        """ u.Vector[length/time, 3]: velocity of this atom; equivalent to
        ``self.momentum/self.mass``
        """
        return (self.momentum / self.mass).defunits()

    @velocity.setter
    def velocity(self, value):
        self.momentum = value * self.mass

    @property
    def num_bonds(self):
        """ int: the number of other atoms this atom is bonded to
        """
        return len(self.bond_graph)
    nbonds = num_bonds

    @property
    def valence(self):
        """ int: the sum of this atom's bond orders
        """
        return sum(v for v in self.bond_graph.itervalues())

    @property
    def symbol(self):
        """ str: elemental symbol
        """
        return data.ELEMENTS.get(self.atnum, '?')
    elem = element = symbol


