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
from ..molecules import MolecularProperties
from ..utils import exports

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
        self._qm_index_set = None

    # TODO: add `qm_atom_indices` to QMMMBase parameters

    def calculate(self, requests):
        self.prep()

        self.mmmol.positions = self.mol.positions
        self._set_qm_positions()

        qmprops = self.qmmol.calculate(requests)
        mmprops = self.mmmol.calculate(requests)

        potential_energy = mmprops.potential_energy+qmprops.potential_energy
        forces = mmprops.forces.copy()
        for iatom, realatom in enumerate(self.qm_atoms):
            forces[realatom.index] = qmprops.forces[iatom]
        for atom in self.qm_link_atoms:
            self._distribute_linkatom_forces(forces, atom)

        properties = MolecularProperties(self.mol,
                                         mmprops=mmprops,
                                         qmprops=qmprops,
                                         potential_energy=potential_energy,
                                         forces=forces)
        if 'wfn' in qmprops:
            properties.wfn = qmprops.wfn

        return properties

    def prep(self):
        if self._prepped:
            return None

        self.params.qm_atom_indices.sort()
        self.qm_atoms = [self.mol.atoms[idx] for idx in self.params.qm_atom_indices]
        self._qm_index_set = set(self.params.qm_atom_indices)

        self.qmmol = self._setup_qm_subsystem()

        self.mmmol = mdt.Molecule(self.mol,
                                  name='%s MM subsystem' % self.mol.name)
        self.mol.ff.copy_to(self.mmmol)
        self._turn_off_qm_forcefield(self.mmmol.ff)
        self.mmmol.set_energy_model(self.params.mm_model)

        self._prepped = True
        return True

    def _setup_qm_subsystem(self):
        raise NotImplemented("%s is an abstract class, use one of its subclasses"
                             % self.__class__.__name__)

    def _turn_off_qm_forcefield(self, ff):
        self._remove_internal_qm_bonds(ff.parmed_obj)
        self._exclude_internal_qm_ljterms(ff.parmed_obj)

    def _exclude_internal_qm_ljterms(self, pmdobj):
        # Turn off QM/QM LJ interactions (must be done AFTER _remove_internal_qm_bonds)
        numqm = len(self.params.qm_atom_indices)
        for i in range(numqm):
            for j in range(i+1, numqm):
                pmdobj.atoms[i].exclude(pmdobj.atoms[j])

    def _remove_internal_qm_bonds(self, pmdobj):
        for i, iatom in enumerate(self.params.qm_atom_indices):
            pmdatom = pmdobj.atoms[iatom]
            allterms = ((pmdatom.bonds, 2), (pmdatom.angles, 3),
                        (pmdatom.dihedrals, 4), (pmdatom.impropers, 4))

            for termlist, numatoms in allterms:
                for term in termlist[:]:  # make a copy so it doesn't change during iteration
                    if self._term_in_qm_system(term, numatoms):
                        term.delete()

    @staticmethod
    def _distribute_linkatom_forces(fullforces, linkatom):
        """ Distribute forces according to the apparently indescribable and unciteable "lever rule"
        """
        # TODO: CHECK THIS!!!!

        mmatom = linkatom.metadata.mmatom
        qmatom = linkatom.metadata.mmpartner
        dfull = mmatom.distance(qmatom)
        d_mm = linkatom.distance(mmatom)

        p = (dfull - d_mm)/dfull
        fullforces[qmatom.index] += p*linkatom.force
        fullforces[mmatom.index] += (1.0-p) * linkatom.force

    def _set_qm_positions(self):
        for qmatom, realatom in zip(self.qmmol.atoms, self.qm_atoms):
            qmatom.position = realatom.position
        mdt.helpers.qmmm.set_link_atom_positions(self.qm_link_atoms)

    def _term_in_qm_system(self, t, numatoms):
        """ Check if an FF term is entirely within the QM subsystem """
        for iatom in range(numatoms):
            attrname = 'atom%i' % (iatom + 1)
            if not getattr(t, attrname).idx in self._qm_index_set:
                return True
        else:
            return False


@exports
class MechanicalEmbeddingQMMM(QMMMEmbeddingBase):
    """
    Handles _non-covalent_ QM/MM with mechanical embedding.

    No electrostatic interactions will be calculated between the QM and MM subsystems.
    No covalent bonds are are allowed between the two susbystems.
    """
    def prep(self):
        if not super(MechanicalEmbeddingQMMM, self).prep():
            return  # was already prepped

        # Set QM partial charges to 0
        self.mmmol.energy_model._prepped = False
        pmdobj = self.mmmol.ff.parmed_obj
        for i, iatom in enumerate(self.params.qm_atom_indices):
            pmdatom = pmdobj.atoms[iatom]
            pmdatom.charge = 0.0

    def _setup_qm_subsystem(self):
        """ QM subsystem for mechanical embedding is the QM atoms + any link atoms
        """
        qm_atoms = [self.mol.atoms[iatom] for iatom in self.params.qm_atom_indices]
        self.qm_link_atoms = mdt.helpers.qmmm.create_link_atoms(self.mol, qm_atoms)
        qmmol = mdt.Molecule(qm_atoms + self.qm_link_atoms,
                             name='%s QM subsystem' % self.mol.name)
        for real_atom, qm_atom in zip(self.qm_atoms, qmmol.atoms):
            qm_atom.metadata.real_atom = real_atom
        qmmol.set_energy_model(self.params.qm_model)
        return qmmol


@exports
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

    def prep(self):
        if not super(ElectrostaticEmbeddingQMMM, self).prep():
            return  # was already prepped

        if not self.params.qm_model.supports_parameter('qm_atom_indices'):
            raise TypeError('Supplied QM model ("%s") does not support QM/MM'
                % self.params.qm_model.__name__)

    def _setup_qm_subsystem(self):
        qmmol = mdt.Molecule(self.mol)
        self.mol.ff.copy_to(qmmol)

        self.qm_link_atoms = mdt.helpers.qmmm.create_link_atoms(self.mol, self.qm_atoms)
        if self.qm_link_atoms:
            raise ValueError('The %s model does not support link atoms' % self.__class__.__name__)

        qmmol.set_energy_model(self.params.qm_model)
        qmmol.energy_model.params.qm_atom_indices = self.params.qm_atom_indices
        return qmmol


