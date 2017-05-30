""" Energy models using OpenBabel's heuristic, highly approximate drug forcefields
"""
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
import numpy as np

import moldesign as mdt
from moldesign import units as u
from moldesign.utils import exports

from .base import EnergyModelBase


# uff disabled because it doesn't pass tests
FFNAMES = {n: n.capitalize() for n in ['ghemical', 'mmff94', 'mmff94s', 'uff']}
UNITNAMES = {'kcal/mol': u.kcalpermol, 'kJ/mol': u.kjpermol}


@exports
class OpenBabelPotential(EnergyModelBase):
    DEFAULT_PROPERTIES = ['potential_energy', 'forces']
    ALL_PROPERTIES = DEFAULT_PROPERTIES
    PARAMETERS = [mdt.parameters.Parameter(name='forcefield',
                                           short_description='Forcefield selection',
                                           type=str,
                                           default='mmff94s',
                                           choices=['ghemical', 'mmff94', 'mmff94s'])]

    def prep(self):
        import openbabel as ob

        if self._prepped: return
        self._pbmol = mdt.interfaces.openbabel.mol_to_pybel(self.mol)
        self._ff = ob.OBForceField.FindForceField(FFNAMES[self.params.forcefield])
        self._ff.Setup(self._pbmol.OBMol)
        self._units = UNITNAMES[self._ff.GetUnit()]

        self._prepped = True

    @mdt.compute.runsremotely(enable=mdt.interfaces.openbabel.force_remote, is_imethod=True)
    def calculate(self, requests):
        import openbabel as ob

        self.prep()
        for atom, pbatom in zip(self.mol.atoms, self._pbmol.atoms):
            pbatom.OBAtom.SetVector(*atom.position.value_in(u.angstrom))
        self._ff.SetCoordinates(self._pbmol.OBMol)
        energy = self._ff.Energy(True)
        self._ff.GetCoordinates(self._pbmol.OBMol)

        # http://forums.openbabel.org/Doing-an-energ-force-calculation-with-openbabel-td1590730.html
        data = ob.toConformerData(self._pbmol.OBMol.GetData(ob.ConformerData))
        force_iterator = data.GetForces()[0]

        forces = np.array([[f.GetX(), f.GetY(), f.GetZ()] for f in force_iterator])
        forces[np.abs(forces) < 1e-30] = 0.0  # prevent underflows

        result = mdt.MolecularProperties(self.mol,
                                         potential_energy=energy*self._units,
                                         forces=forces*self._units/u.angstrom)
        return result


    # if necessary, run the entire minimization remotely for speed
    @mdt.compute.runsremotely(enable=mdt.interfaces.openbabel.force_remote, is_imethod=True)
    def minimize(self, **kwargs):
        super().minimize(**kwargs)
