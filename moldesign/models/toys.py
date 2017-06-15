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
from .. import units as u
from ..utils import exports
from .base import EnergyModelBase


@exports
class Spring(EnergyModelBase):
    """Two atoms attached by a spring"""

    DEFAULT_PROPERTIES = ['potential_energy',
                          'force']
    ALL_PROPERTIES = DEFAULT_PROPERTIES

    PARAMETERS = [
        mdt.parameters.Parameter(
            'k', 'Spring constant',
            type=u.eV/u.angstrom**2),
        mdt.parameters.Parameter(
            'd0', 'Equilibrium distance',
            type=u.angstrom)]

    def prep(self):
        assert self.mol.natoms == 2
        self._prepped = True

    def calculate(self, requests):
        dvec = self.mol.atoms[1].position - self.mol.atoms[0].position
        dist = np.sqrt(dvec.dot(dvec))
        pe = 0.5 * self.params.k * (dist - self.params.d0)**2
        f = self.params.k * dvec * (dist - self.params.d0) / dist
        forcearray = u.array([f, -f])
        return {'potential_energy': pe,
                'forces': forcearray}


@exports
class HarmonicOscillator(EnergyModelBase):
    """ Applies a harmonic potential (centered at 0) to the x-component of every atom
    """
    PARAMETERS = [
        mdt.parameters.Parameter(
            'k', 'Spring constant',
            type=u.eV/u.angstrom**2)]

    def prep(self):
        if self.params.k.dimensionality != (u.eV/u.angstrom**2).dimensionality:
            raise u.DimensionalityError("Spring constant must have dimensions of energy/length^2")
        return True

    def calculate(self, requests):
        self.prep()
        energy = 0.5 * self.params.k * np.sum(self.mol.positions[:, 0]**2)
        forces = np.zeros((self.mol.num_atoms, 3)) * u.default.force
        forces[:, 0] = - self.params.k * self.mol.positions[:, 0]
        return dict(potential_energy=energy, forces=forces)
