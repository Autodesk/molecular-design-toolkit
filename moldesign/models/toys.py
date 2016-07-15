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

from moldesign import units as u
from .base import EnergyModelBase


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class Spring(EnergyModelBase):
    """Two atoms attached by a spring"""
    def __init__(self, d0, k):
        self.d0 = d0
        self.k = k
        super(Spring, self).__init__()

    def prep(self):
        assert self.mol.natoms == 2
        self._prepped = True

    def calculate(self, requests):
        dvec = self.mol.atoms[1].position - self.mol.atoms[0].position
        dist = np.sqrt(dvec.dot(dvec))
        pe = 0.5 * self.k * (dist - self.d0)**2
        f = self.k * dvec * (dist - self.d0) / dist
        forcearray = u.array([f, -f])
        return {'potential_energy': pe,
                'forces': forcearray}


@exports
class HarmonicOscillator(EnergyModelBase):
    """ Applies a harmonic potential (centered at 0) to the x-component of every atom
    """
    def __init__(self, k):
        self.k = k
        assert k.dimensionality == {'[mass]': 1,
                                    '[time]': -2}
        super(EnergyModelBase,self).__init__()

    def calculate(self, requests):
        energy = 0.5 * self.k * np.sum(self.mol.positions[:, 0]**2)
        forces = np.zeros((self.mol.num_atoms, 3)) * u.hartree / u.bohr
        forces[:, 0] = - self.k * self.mol.positions[:, 0]
        return dict(potential_energy=energy, forces=forces)