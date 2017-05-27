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
import json

import pyccc

import moldesign as mdt
from moldesign import units as u
from moldesign.utils import exports

from .base import QMBase

IMAGE = 'nwchem'


@exports
class NWChemQM(QMBase):
    """ Interface with NWChem package (QM only)

    Note:
        This is the first interface based on our new wrapping strategy. This is slightly hacked,
        but has the potential to become very general; very few things here are NWChem-specific
    """

    DEFAULT_PROPERTIES = ['potential_energy']

    def prep(self):
        parameters = self.params.copy()

        parameters['constraints'] = []
        for constraint in self.mol.constraints:
            raise NotImplementedError()

        parameters['charge'] = self.mol.charge.value_in(u.q_e)
        self._jobparams = parameters
        self._prepped = True

    def calculate(self, requests=None):
        if requests is None: requests = self.DEFAULT_PROPERTIES
        job = self._make_calculation_job(requests)

        return mdt.compute.run_job(job, _return_result=True)

    def _make_calculation_job(self, requests=None):
        self.prep()
        parameters = self._jobparams.copy()
        parameters['runType'] = 'singlePoint'
        parameters['properties'] = list(requests)
        if self.mol.constraints:
            self.write_constraints(parameters)
        job = pyccc.Job(  # image=mdt.compute.get_image_path(IMAGE),
                image=IMAGE,
                command='run.py && getresults.py',
                inputs={'input.xyz': self.mol.write(format='xyz'),
                        'params.json': json.dumps(parameters)},
                when_finished=self.finish,
                name='nwchem/%s'%self.mol.name)
        return job

    def write_constraints(self, parameters):
        parameters['constraints'] = []
        for constraint in self.mol.constraints:

            if not constraint.satisfied():  # TODO: factor this out into NWChem subclass
                raise ValueError('Constraints must be satisfied before passing to NWChem %s'
                                 %constraint)

            cjson = {'type': constraint['desc'],
                     'value': constraint['value'].to_json()}
            if cjson['type'] == 'position':
                cjson['atomIdx'] = constraint.atom.index
            elif cjson['type'] == 'distance':
                cjson['atomIdx1'] = constraint.a1.index
                cjson['atomIdx2'] = constraint.a2.index


    def finish(self, job):
        results = json.loads(job.get_output('results.json').read())
        return self._process_results(results)

    def _process_results(self, results):
        assert len(results['states']) == 1
        jsonprops = results['states'][0]['calculated']
        result = mdt.MolecularProperties(self.mol,
                                         **self._get_properties(jsonprops))
        return result

    @staticmethod
    def _get_properties(jsonprops):
        props = {}
        for name, property in jsonprops.items():
            if isinstance(property, dict) and len(property) == 2 and \
                            'units' in property and 'value' in property:
                props[name] = property['value'] * u.ureg(property['units'])
            else:
                props[name] = property

        return props

    def minimize(self, nsteps=None):
        job = self._make_minimization_job(nsteps)

        return mdt.compute.run_job(job, _return_result=True)

    def _make_minimization_job(self, nsteps):
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
        return job

    def finish_min(self, job):
        # TODO: parse more data than just the final minimization state
        properties = self.finish(job)
        traj = mdt.Trajectory(self.mol)
        traj.new_frame()
        self.mol.positions = properties.positions
        self.mol.properties = properties
        traj.new_frame()
        return traj


