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

from ..compute import runsremotely
from ..interfaces import psi4_interface
from .base import QMBase

THEORIES = {'rhf':'scf'}

BASES = {}


class Psi4Potential(QMBase):
    def prep(self):
        return True

    @runsremotely(enable=psi4_interface.force_remote)
    def calculate(self, requests):
        import psi4
        geom = psi4_interface.mdt_to_psi4(self.mol)
        psi4.set_options({'basis': BASES.get(self.params.basis, self.params.basis)})
        return {'potential_energy': self._get_runmethod(requests)(self._gettheorystring())}

    @staticmethod
    def _get_runmethod(requests):
        import psi4
        # TODO: actually deal with requests
        return psi4.energy

    def _gettheorystring(self):
        return THEORIES.get(self.params.theory, self.params.theory)



