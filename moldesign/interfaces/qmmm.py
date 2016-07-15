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
from moldesign.methods import basemethods


class NonbondedQMMM(basemethods.QMMMBase):
    """
    A simple QM/MM prescription that doesn't treat bonds between the QM and MM systems.
    The QM and MM subsystems are coupled via Lennard-Jones and coulombic interactions:
    $$H = H_{QM} + H_{MM} + H{inter},$$ where
    $$H_{inter} = H_{LJ} + H_{coulomb}$$.
    In practical terms, this means that the QM calculation must include the MM atoms as point charges,
    while the MM calculation must include the QM atoms as unbonded, neutral lennard-jones spheres.
    """

    def prep(self):
        for atom in self.params.qm_atoms: atom.props.subsystem = 'qm'
        for atom in self.params.mm_atoms: atom.props.subsystem = 'mm'

        self._build_mm_system()
        self._build_qm_system()

        self._prepped = True


    def _build_qm_system(self):
        # Set up the QM system
        self.qmmol = mdt.Molecule(self.mol)
        for atom, subatom in zip(self.mol.atoms, self.qmmol.atoms):
            subatom.molecule_atom = atom
            atom.props.qm_child = subatom

        # Make the MM atoms point charges
        point_charges = {atom: atom.props.mm_child.ff.charge for atom in self.qmmol.atoms
                         if atom.molecule_atom.params.subsystem == 'mm'}
        self.qm_model.params.point_charges = point_charges
        self.qmmol.set_energy_model(self.qm_model)
        self._subatoms = self.qmmol.atoms + self.mmol.atoms


    def _build_mm_system(self):
        # Set up the MM subsystem - includes all atoms,
        # but we remove all FF terms except LJ
        self.mmmol = mdt.Molecule(self.mol)
        for atom, subatom in zip(self.mol.atoms, self.mmmol.atoms):
            subatom.molecule_atom = atom
            atom.props.mm_child = subatom
        self.mmmol.set_energy_model(self.mm_model)
        # Remove charges and bonds for QM atoms
        ff = self.mm_model.params.ff
        forcefield = ff.get_parameters(self.mmmol)

        def prune(termlist):  # remove intra-qm region forces
            # TODO: what about LJ forces???
            newterms = []
            for term in termlist:
                subsystems = set([atom.molecule_atom.props.subsystem for atom in termlist.atom])
                subsystems = list(subsystems)
                if len(subsystems) == 1 and subsystems[0] == 'mm':
                    newterms.append(subsystems)
                elif len(subsystems) > 1:
                    raise ValueError('Bond crosses boundary: %s' % term)

        forcefield.bonds = prune(forcefield.bonds)
        forcefield.angles = prune(forcefield.angles)
        forcefield.dihedrals = prune(forcefield.dihedrals)
        forcefield.impropers = prune(forcefield.impropers)

    def calculate(self, requests=None):
        for atom in self._subatoms:
            atom.position = atom.molecule_atom.position

        qmprops = self.qm_model.calculate()
        mmprops = self.mm_model.calculate()

