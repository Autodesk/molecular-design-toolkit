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

import operator as op

import moldesign as mdt
from moldesign.parameters import Parameter, WhenParam

from . import ForceField


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class GAFF(ForceField):
    """ Model the energy using the GAFF forcefield

    This is implemented as a special case of the ForceField energy model; it automates small
    parameterization process
    """
    # TODO: mechanism to store partial charges so they don't need to be constantly recomputed

    PARAMETERS = [Parameter('partial_charges',
                            'Partial charge model',
                            type=str,
                            default='am1-bcc',
                            choices=['am1-bcc', 'gasteiger', 'esp']),
                  Parameter('gaff_version',
                            'GAFF version',
                            type=str,
                            choices='gaff gaff2'.split(),
                            default='gaff2')
                  ] + ForceField.PARAMETERS

    def prep(self, force=False):
        self._parameterize()
        return super(GAFF, self).prep()

    def calculate(self, requests=None):
        if not self._prepped:
            self._parameterize()
        return super(GAFF, self).calculate(requests=requests)

    def _parameterize(self):
        if not self.mol.ff:
            mdt.parameterize(self.mol,
                             charges=self.params.partial_charges,
                             ffname=self.params.gaff_version)

