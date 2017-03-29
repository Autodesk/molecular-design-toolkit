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
from itertools import chain as iterchain

from moldesign import utils
from moldesign import units as u
from .atomcollections import AtomListOperationMixin


class BaseMasterList(utils.AutoIndexList):
    """ Base class for a list of objects (atoms, residues, etc.) in a molecule
    """
    def __init__(self, objs, mol):
        super(BaseMasterList, self).__init__(objs)
        self._mol = mol

    def _error_if_object_owned(self, obj):
        if obj.molecule is self._mol:
            raise ValueError("%s cannot appear twice in this molecule" % obj)
        elif obj.molecule is not None:
            raise ValueError("%s is already a member of %s" % (obj, obj.molecule))

    def append(self, obj):
        """ Adds obj to this molecule. Equivalent to ``obj.molecule = mol``

        Args:
            obj (object): object to add

        Note:
            It's much more efficient to use the ``extend`` method for lists of objects

        Raises:
            ValueError: If obj is already part of a different molecule, or is part of
               a larger substructure (such as a chain or residue)
        """
        self._error_if_object_owned(obj)
        super(BaseMasterList, self).append(obj)

    def extend(self, objs):
        """ Adds a list of objects to this molecule. More efficient than adding them individually.

        If the object is already assigned to a residue, the residue will be assigned to this molecule
        as well.

        Args:
            obj (List[object]): objects to add

        Raises:
            ValueError: If obj is already part of a different molecule
        """
        for o in objs:
            self._error_if_object_owned(o)
        super(BaseMasterList, self).extend(objs)

    def insert(self, index, obj):
        self._error_if_object_owned(obj)
        super(BaseMasterList, self).insert(index, obj)
    insert.__doc__ = append.__doc__

    def remove(self, obj):
        raise NotImplementedError()  # must be implemented by subclass

    def _remove_from_list(self, obj):
        super(BaseMasterList, self).remove(obj)

    def pop(self, index=-1):
        obj = self[index]
        self.remove(obj)
        return obj
    pop.__doc__ = remove.__doc__


class AtomMasterList(BaseMasterList, AtomListOperationMixin):
    """ The master list of atoms in a molecule.
    """
    def _error_if_obj_owned(self, atom):
        if atom.residue is not None and atom.residue.molecule is not self._mol:
            raise ValueError("Cannot individually assign %s to %s, because it is part of %s" %
                             (atom, self._mol, atom.residue))
        super(AtomMasterList, self)._error_if_object_owned(atom)

    #@utils.doc_inherit
    def append(self, atom):
        self.extend([atom])

    def extend(self, atoms):
        """
        Note:
            Internally, this is THE place where new atoms are added to a molecular structure
        """
        # Do error checking immediately so that we don't get to an inconsistent state.
        # Any errors after this point will probably screw up the molecule
        numinitial = len(self)
        self._error_if_cant_add_atoms(atoms)
        for o in atoms:
            self._error_if_object_owned(o)

        # Resize the physics arrays
        self._mol.ndims = 3 * len(self) + len(atoms)
        self._mol.masses = u.array([atom.mass for atom in iterchain(self, atoms)]).defunits()
        self._mol.dim_masses = u.broadcast_to(self._mol.masses, (3, len(self)+len(atoms))).T
        self._mol._dof = None
        self._mol._positions = u.array([atom.position for atom in iterchain(self, atoms)])
        self._mol._momenta = u.array([atom.momentum for atom in iterchain(self, atoms)])

        # assign default residues and chains where necessary
        for iatom, atom in enumerate(atoms):
            utils.AutoIndexList.append(self, atom)

            if atom.residue is None:
                self._mol._defresidue.add(atom, _addatoms=False)

            if atom.residue.chain is None:
                self._mol._defchain.add(atom.residue, _addatoms=False)

            if atom.chain.molecule is not self._mol:
                assert atom.chain.molecule is None
                self._mol.chains.add(atom.chain, _addatoms=False)

            self._mol.bonds._addatom(atom)
            assert atom.residue.molecule is self._mol
            assert atom.residue in self._mol.residues
            assert atom.chain.molecule is self._mol
            assert atom.chain in self._mol.chains
            assert atom in self._mol.atoms

        # Mark the atoms as owned
        for iatom, atom in enumerate(atoms):
            atom._delegate_state_to_molecule(self._mol)

    def _error_if_cant_add_atoms(self, atoms):
        atoms_to_add = set(atoms)
        structure_to_add = set()
        atoms_in_structure = set()
        for atom in atoms:
            for obj in (atom.residue, atom.chain):
                if obj is not None and obj not in structure_to_add:
                    structure_to_add.add(obj)
                    atoms_in_structure.update(obj.atoms)

            for nbr in atom.bonds.atoms:
                if nbr.molecule is not None and nbr.molecule is not self._mol:
                    raise ValueError(
                            "Can't add atoms to %s: %s is bonded to a different molecule ('%s')" %
                            (self._mol, atom, nbr.molecule))

        if atoms_in_structure > atoms_to_add:
            raise ValueError("Can't add atoms - they're part of a larger structure")

    @utils.doc_inherit
    def pop(self, index=-1):
        atom = self[index]
        self.remove(atom)
        return atom

    def remove(self, atom):
        """ Removes ``atom`` from all larger structures (molecule, chain, and residue).

        Equivalent to ``obj.residue = None``

        Args:
            obj (object): object to remove

        Raises:
            ValueError: If obj is not part of this
        """
        atom._recover_state_from_molecule()
        atom.residue._remove(atom)  # this removes the atom

    def _remove_from_list_and_bonds(self, atom):
        self._mol.bonds._removeatom(atom)
        self._remove_from_list(atom)

    def _append_and_update_bonds(self, atom):
        self._mol.bonds._addatom(atom)
        utils.AutoIndexList.append(self, atom)

    def _insert_and_update_bonds(self, index, atom):
        self._mol.bonds._addatom(atom)
        utils.AutoIndexList.insert(self, index, atom)

    def _extend_and_update_bonds(self, atoms):
        for atom in atoms:
            self._mol.bonds._addatom(atom)
        utils.AutoIndexList.extend(self, atoms)


class ResidueMasterList(BaseMasterList):
    """ The master list of residues in a molecule.
    """
    def _error_if_obj_owned(self, res):
        if res.chain is not None and res.molecule is not self._mol:
            raise ValueError("Cannot individually assign %s to %s, because it is part of %s" %
                             (res, self._mol, res.residue))
        super(ResidueMasterList, self)._error_if_object_owned(res)

    @utils.doc_inherit
    def append(self, res):
        self.extend([res])

    def pop(self, index=-1):
        res = super(ResidueMasterList, self).pop(index)
        res._chain = None
        return res

    def extend(self, residues):
        atoms = []
        for r in residues:
            atoms.extend(r.atoms)
        self._mol.atoms.extend(atoms)

    def remove(self, residue):
        assert residue is self[residue.index]
        residue.chain._remove(residue)
