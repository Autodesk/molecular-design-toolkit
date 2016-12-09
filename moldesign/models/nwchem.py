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

def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class NWChemQM(QMBase):
    """ Interface with NWChem package (QM only)

    Note:
        This is the first interface based on our new wrapping strategy. This is slightly hacked,
        but has the potential to become very general; very few things here are NWChem-specific
    """

    def prep(self):
        parameters = self.params.copy()

        parameters['constraints'] = []
        for constraint in self.mol.constraints:
            raise NotImplementedError()

        parameters['charge'] = self.mol.charge.value_in(u.q_e)
        self._jobparams = parameters
        self._prepped = True

    def calculate(self, requests=None):
        self.prep()

        parameters = self._jobparams.copy()
        parameters['runType'] = 'singlePoint'
        parameters['properties'] = list(requests)

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

    def minimize(self, nsteps=None):
        self.prep()

        parameters = self._jobparams.copy()

        parameters['runType'] = 'minimization'
        if nsteps is not None:
            parameters['minimization_steps'] = 100

        job = pyccc.Job(  # image=mdt.compute.get_image_path(IMAGE),
                image=IMAGE,
                command='run.py && getresults.py',
                inputs={'input.xyz': self.mol.write(format='xyz'),
                        'params.json': json.dumps(parameters)},
                when_finished=self.finish_min,
                name='nwchem/%s' % self.mol.name)

        return mdt.compute.run_job(job, _return_result=True)

    def finish_min(self, job):
        # TODO: do a better job here
        properties = self.finish(job)
        traj = mdt.Trajectory(self.mol)
        traj.new_frame()
        self.mol.positions = properties.positions
        self.mol.properties = properties
        traj.new_frame()
        return traj


