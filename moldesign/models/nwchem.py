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
import json

import pyccc

import moldesign as mdt
from moldesign import units as u

from .base import QMBase

IMAGE = 'nwchem'

class NWChemQM(QMBase):
    """ Interface with NWChem package (QM only)

    Note:
        This is the first interface based on our new wrapping strategy. This is slightly hacked,
        but has the potential to become very general; very few things here are NWChem-specific
    """

    def prep(self):
        self._prepped = True

    def calculate(self, requests=None):
        self.prep()

        parameters = self.params.copy()
        if requests is not None:
            parameters['properties'] = list(requests)

        parameters['constraints'] = []
        for constraint in self.mol.constraints:
            raise NotImplementedError()

        parameters['charge'] = self.mol.charge.value_in(u.q_e)
        parameters['runType'] = 'singlePoint'

        job = pyccc.Job(  # image=mdt.compute.get_image_path(IMAGE),
                image=IMAGE,
                command='run.py && getresults.py',
                inputs={'input.xyz': self.mol.write(format='xyz'),
                        'params.json': json.dumps(parameters)},
                when_finished=self.finish,
                name='nwchem/%s' % self.mol.name)

        return mdt.compute.run_job(job, _return_result=True)

    def finish(self, job):
        results = json.loads(job.get_output('results.json').read())
        assert len(results['states']) == 1

        jsonprops = results['states'][0]['calculated']
        result = mdt.MolecularProperties(self.mol,
                                         **self._get_properties(jsonprops))
        return result

    @staticmethod
    def _get_properties(jsonprops):
        props = {}
        for name, property in jsonprops.iteritems():
            if isinstance(property, dict) and len(property) == 2 and \
                            'units' in property and 'value' in property:
                props[name] = property['value'] * u.ureg(property['units'])
            else:
                props[name] = property

        return props

