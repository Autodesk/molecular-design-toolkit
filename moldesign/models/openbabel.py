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
""" Energy models using OpenBabel's heuristic, highly approximate drug forcefields
"""
from __future__ import absolute_import

import openbabel as ob

import moldesign as mdt
from moldesign import units as u
from .base import EnergyModelBase


class OpenBabelPotential(EnergyModelBase):
    DEFAULT_PROPERTIES = ['potential_energy', 'forces']
    ALL_PROPERTIES = DEFAULT_PROPERTIES
    PARAMETERS = [mdt.parameters.Parameter(name='forcefield',
                                           short_description='Forcefield selection',
                                           type=str,
                                           default='mmff94s',
                                           choices=['ghemical', 'mmff94', 'mmff94s', 'uff'])]
    FORCE_UNITS = u.hartree / u.bohr

    FFNAMES = {n: n.capitalize() for n in ['ghemical', 'mmff94', 'mmff94s', 'uff']}

    def prep(self):
        if self._prepped:
            return
        self._pbmol = mdt.interfaces.openbabel.mol_to_pybel(self.mol)
        self._ff = ob.OBForceField.FindForceField(self.FFNAMES[self.params.forcefield])
        self._ff.Setup(self._pbmol.OBMol)

        assert self._ff.GetUnit() == 'kcal/mol'
        self._prepped = True

    def calculate(self, requests):
        self.prep()
        for atom, pbatom in zip(self.mol.atoms, self._pbmol.atoms):
            pbatom.OBAtom.SetVector(*atom.position.value_in(u.angstrom))
        self._ff.UpdateCoordinates(self._pbmol.OBMol)

        energy = self._ff.Energy(True)
        data = ob.toConformerData(self._pbmol.OBMol.GetData(4))
        forces = []

        for iatom, f in enumerate(data.GetForces()[0]):
            forces.append([f.GetX(), f.GetY(), f.GetZ()])

        result = mdt.MolecularProperties(self.mol,
                                         potential_energy=energy*u.kcalpermol,
                                         forces=forces*u.kcalpermol/u.angstrom)

        return result






