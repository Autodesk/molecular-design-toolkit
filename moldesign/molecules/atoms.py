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
from .notebook_display import AtomNotebookMixin


class AtomPropertyMixin(object):
    """ Functions accessing computed atomic properties.

    Note:
        This is a mixin class designed only to be mixed into the :class:`Atom` class. Routines
        are separated are here for code organization only - they could be included in the main
        Atom class without changing any functionality
    """
    distance = utils.Alias('_container.distance')
    atoms_within = utils.Alias('_container.atoms_within')
    residues_within = utils.Alias('_container.residues_within')

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

    @property
    def ff(self):
        """ moldesign.utils.DotDict: This atom's force field parameters, if available (``None``
        otherwise)
        """
        if self.molecule.ff is None:
            return None
        else:
            return self.molecule.ff.get_atom_terms(self)

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


def _on_name_change(atom, oldname, newname):
    if getattr(atom, 'residue', None):
        atom.residue._remove(atom)
        atom.residue._add(atom)


@toplevel
class Atom(AtomPropertyMixin, AtomNotebookMixin):
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
            - :class:`AtomNotebookMixin`
    """
    x, y, z = (AtomCoordinate('position', i) for i in xrange(3))
    vx, vy, vz = (AtomCoordinate('velocity', i) for i in xrange(3))
    px, py, pz = (AtomCoordinate('momentum', i) for i in xrange(3))
    fx, fy, fz = (AtomCoordinate('force', i) for i in xrange(3))
    position = AtomArray('_position', 'positions')
    momentum = AtomArray('_momentum', 'momenta')

    name = utils.EventfulAttr('_name', _on_name_change)

    atomic_number = utils.Synonym('atnum')

    #################################################################
    # Methods for BUILDING the atom and indexing it in a molecule
    def __init__(self, name=None, atnum=None, mass=None, residue=None,
                 formal_charge=None, pdbname=None, pdbindex=None, element=None,
                 metadata=None):

        # Allow user to instantiate an atom as Atom(6) or Atom('C')
        if atnum is None and element is None:
            if isinstance(name, int):
                atnum = name
                name = None
            else:
                element = name

        if element:
            self.atnum = data.ATOMIC_NUMBERS[element]
        else:
            self.atnum = atnum

        self.name = utils.if_not_none(name, self.elem)
        self.pdbname = utils.if_not_none(pdbname, self.name)
        self.pdbindex = pdbindex

        if mass is None: self.mass = data.ATOMIC_MASSES[self.atnum]
        else: self.mass = mass

        self.formal_charge = utils.if_not_none(formal_charge, 0.0 * u.q_e)
        self._residue = residue
        self._index = None
        self._position = np.zeros(3) * u.default.length
        self._momentum = np.zeros(3) * (u.default.length * u.default.mass/u.default.time)
        self._bond_graph = {}
        self.metadata = mdt.utils.DotDict()
        if metadata:
            self.metadata.update(metadata)

    @property
    def index(self):
        return self._index

    @property
    def molecule(self):
        if self.residue and self.residue.chain:
            return self.residue.chain.molecule
        else:
            return None

    @molecule.setter
    def molecule(self, mol):
        if mol is self.molecule:
            return
        if self.molecule is not None:
            self.molecule.atoms.remove(self)
        if mol is not None:
            mol.atoms.append(self)
            assert self.molecule is mol
            assert self.chain is mol._defchain
            assert self.residue is mol._defresidue

    @property
    def residue(self):
        return self._residue

    @residue.setter
    def residue(self, res):
        if res is self._residue:
            return
        if self.molecule is not None:  # pop it out of the molecule now
            self._molecule.atoms.remove(self)
        assert self.molecule is None
        if res is not None:
            res.add(self)

    @property
    def chain(self):
        if self.residue:
            return self.residue.chain
        else:
            return None

    def _add_to_structure(self, struc, attr):
        if getattr(self, attr) is struc:
            return
        if getattr(self, attr) is not None:
            self.molecule.remove(self)
        if struc is not None:
            getattr(self, attr)._add(self)

    def __str__(self):
        desc = '%s %s (elem %s)' % (self.__class__.__name__, self.name, self.elem)
        molstring = ''
        if self.molecule:
            molstring = ', index %s' % self.index
            if self.molecule.is_biomolecule:
                molstring += ' (res %s chain %s)' % (self.residue.name, self.chain.name)
        return '%s%s' % (desc, molstring)

    def _shortstr(self):
        """ A shorter string representation for easier-to-read lists of atoms
        """
        fields = [self.name]
        if self.molecule:
            fields.append('#%d' % self.index)
            if self.molecule.is_biomolecule:
                fields.append('in %s.%s' % (self.chain.name, self.residue.name))
        return ' '.join(fields)

    def __repr__(self):
        try:
            if self.molecule:
                return '<%s in molecule %s>' % (self, self.molecule)
            else:
                return '<%s>' % self
        except (KeyError, AttributeError):
            return '<%s at %s (exception in __repr__)>' % (self.__class__.__name__, id(self))

    @utils.args_from(AtomContainer.copy_atoms)
    def copy(self, *args, **kwargs):
        """ Copy an atom (delegate to AtomContainer)
        """
        return self._container.copy_atoms(*args, **kwargs)[0]

    def __getstate__(self):
        """Helper for pickling"""
        state = self.__dict__.copy()
        if self.molecule is not None:  # then these don't belong to the atom anymore
            state['_bond_graph'] = None
            state['_position'] = self.position
            state['_momentum'] = self.momentum
        return state

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
    def bonded_atoms(self):
        """ List[moldesign.Atom]: a list of the atoms this atom is bonded to
        """
        return [bond.partner(self) for bond in self.bonds]

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


