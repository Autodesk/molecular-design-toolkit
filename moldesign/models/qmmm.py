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

import moldesign as mdt
from moldesign import units as u
from moldesign.molecules import MolecularProperties

from .base import QMMMBase


class QMMMEmbeddingBase(QMMMBase):
    """ Abstract class for standard QM/MM embedding models.

    To use any of this classes' subclasses, the MM models must support the ability to
    calculate the internal energies and interaction energies between subsystems,
    using the ``calculation_groups`` parameter.
    """

    def __init__(self, *args, **kwargs):
        super(QMMMEmbeddingBase, self).__init__(*args, **kwargs)
        self.qmmol = None
        self.mmmol = None
        self.qm_atoms = None
        self.qm_link_atoms = None

    # TODO: add `qm_atom_indices` to QMMMBase parameters
    # TODO: add `subsystem_atom_indices` to MMBase parameters

    def calculate(self, requests):
        self.mmmol.positions = self.mol.positions
        self._set_qm_positions()

        qmprops = self.qmmol.calculate(requests)
        mmprops = self.mmmol.calculate(requests)

        properties = MolecularProperties(self.mol,
                                         mmprops=mmprops,
                                         qmprops=qmprops)

        for additive_prop in ('potential_energy', 'forces'):
            properties['additive_prop'] = (qmprops[additive_prop] +
                                           mmprops[additive_prop] -
                                           mmprops.qm_system[additive_prop])

    def prep(self):
        if self._prepped:
            return None

        self.params.qm_atom_indices.sort()
        self.qm_atoms = [self.mol.atoms[idx] for idx in self.params.qm_atom_indices]
        self._qm_index_set = set(self.params.qm_atom_indices)

        self.qmmol = self._setup_qm_subsystem()
        self.mmmol = mdt.Molecule(self.mol,
                                  name='%s MM subsystem' % self.mol.name)
        # Disabled due to https://github.com/ParmEd/ParmEd/issues/839
        # self.mol.ff.copy_to(self.mmmol)
        self.mmmol.ff = mdt.forcefields.ForceField(self.mmmol, self.mol.ff.sourcedata)

        self.mmmol.set_energy_model(self.params.mm_model, **self.params)
        self._separate_qm_and_mm_systems(self.mmmol.energy_model)

        self._prepped = True
        return True

    def _setup_qm_subsystem(self):
        raise NotImplemented("%s is an abstract class, use one of its subclasses"
                             % self.__class__.__name__)

    def _set_qm_positions(self):
        raise NotImplemented("%s is an abstract class, use one of its subclasses"
                             % self.__class__.__name__)

    def _separate_qm_and_mm_systems(self, model):
        calculation_groups = model.params.setdefault('calculation_groups', {})
        assert 'qm_system' not in calculation_groups
        assert 'mm_system' not in calculation_groups

        calculation_groups['qm_system'] = sorted(self.params.qm_atom_indices)
        calculation_groups['mm_system'] = [i for i in xrange(self.mol.num_atoms)
                                           if i not in self._qm_index_set]


class MechanicalEmbeddingQMMM(QMMMEmbeddingBase):
    """
    Handles _non-covalent_ QM/MM with mechanical embedding.

    No electrostatic interactions will be calculated between the QM and MM subsystems.
    No covalent bonds are are allowed between the two susbystems.
    """

    def prep(self):
        if not super(MechanicalEmbeddingQMMM, self).prep():
            return  # was already prepped

        # Set forcefield partial charges for the QM subsystem to 0
        self.mmmol.energy_model._prepped = False
        for atom in self.qm_atoms:
            atom.ff.partial_charge = 0.0 * u.q_e

    def _set_qm_positions(self):
        for atom in self.qmmol.atoms:
            atom.position = atom.metadata.real_atom.position
            mdt.helpers.qmmm.set_link_atom_positions(self.qm_link_atoms)

    def _setup_qm_subsystem(self):
        """ QM subsystem for mechanical embedding is the QM atoms + any link atoms
        """
        qmatoms = [self.mol.atoms[iatom] for iatom in self.params.qm_atom_indices]
        self.qm_link_atoms = mdt.helpers.qmmm.create_link_atoms(self.mol, qmatoms)
        qmmol = mdt.Molecule(qmatoms + self.qm_link_atoms,
                             name='%s QM subsystem' % self.mol.name)
        for real_atom, qm_atom in zip(self.qm_atoms, qmmol.atoms):
            qm_atom.metadata.real_atom = real_atom
        qmmol.set_energy_model(self.params.qm_model, **self.params)
        return qmmol


class ElectrostaticEmbeddingQMMM(QMMMEmbeddingBase):
    """ Handles _non-covalent_ QM/MM with electrostaic embedding.
    No bonds allowed across the QM/MM boundaries.

    To support this calculation type, the QM model must support the ability to denote
    a subset of atoms as the "QM" atoms, using the ``qm_atom_indices`` parameter. 

    To support this calculation type, the QM model must support the ability to denote
    a subset of atoms as the "QM" atoms, using the ``qm_atom_indices`` parameter. 
    The MM models must support the ability to turn of _internal_ interactions for 
    a certain subset of the system, using the ``no_internal_calculations`` parameter.

    """

    PARAMETERS = QMMMBase.PARAMETERS + ['qm_model', 'mm_model']

    def prep(self):
        if not super(ElectrostaticEmbeddingQMMM, self).prep():
            return  # was already prepped

        if 'qm_atom_indices' not in self.params.qm_model.PARAMETERS:
            raise TypeError('Supplied QM model ("%s") does not support QM/MM'
                % self.params.qm_model.__name__)



    def _setup_qm_subsystem(self):
        qmmol = self.mol.copy()
        qmmol.set_energy_model(self.params.qm_model, 
                               qm_atom_indices=list(self.params.qm_atom_indices),
                               **self.params)
        return qmmol






