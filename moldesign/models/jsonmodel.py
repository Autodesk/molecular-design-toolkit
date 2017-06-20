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

from .base import EnergyModelBase


class JsonModelBase(EnergyModelBase):
    """ Abstract energy model interface using JSON inputs and outputs
    """
    IMAGE = None  # base name of the docker image
    MODELNAME = 'abstract model'

    RUNNER = 'run.py'
    PARSER = 'getresults.py'

    def prep(self):
        parameters = self.params.copy()

        parameters['constraints'] = []
        for constraint in self.mol.constraints:
            self._handle_constraint(constraint)

        parameters['charge'] = self.mol.charge.value_in(u.q_e)
        self._jobparams = parameters
        self._prepped = True

    def calculate(self, requests=None):
        if requests is None:
            requests = self.DEFAULT_PROPERTIES
        job = self._make_calculation_job(requests)

        return mdt.compute.run_job(job, _return_result=True)

    def _handle_constraint(self, constraint):
        raise NotImplementedError()

    def _make_calculation_job(self, requests=None):
        params, inputfiles = self._prep_calculation(requests)
        inputfiles['params.json'] = mdt.utils.json_dumps(dict(params))
        job = pyccc.Job(image=mdt.compute.get_image_path(self.IMAGE),
                        command='%s && %s' % (self.RUNNER, self.PARSER),
                        inputs=inputfiles,
                        when_finished=self.finish,
                        name='%s/%s' % (self.MODELNAME, self.mol.name))
        return job

    def _prep_calculation(self, requests):
        self.prep()
        parameters = self._jobparams.copy()
        parameters['runType'] = 'singlePoint'
        parameters['properties'] = list(requests)
        if self.mol.constraints:
            self.write_constraints(parameters)
        inputfiles = self._get_inputfiles()
        return parameters, inputfiles

    def finish(self, job):
        results = json.loads(job.get_output('results.json').read())
        return self._process_results(results)

    def _process_results(self, results):
        assert len(results['states']) == 1
        jsonprops = results['states'][0]['calculated']
        if 'orbitals' in jsonprops:
            wfn = self._make_wfn(results['states'][0])
        else:
            wfn = None

        result = mdt.MolecularProperties(self.mol,
                                         **self._json_to_quantities(jsonprops))
        if wfn:
            result['wfn'] = wfn
        return result

    def _make_wfn(self, state):
        from moldesign import orbitals

        try:
            basis_fns = state['calculated']['method']['aobasis']
        except KeyError:
            basis_set = None
        else:
            bfs = [orbitals.AtomicBasisFunction(**bdata) for bdata in basis_fns]
            basis_set = orbitals.BasisSet(self.mol,
                                          orbitals=bfs,
                                          name=self.params.basis)
        wfn = orbitals.ElectronicWfn(self.mol,
                                     self.mol.num_electrons,
                                     aobasis=basis_set)

        for setname, orbdata in state['calculated']['orbitals'].items():
            orbs = []
            for iorb in range(len(orbdata['coefficients'])):
                orbs.append(orbitals.Orbital(orbdata['coefficients'][iorb]))
                if 'occupations' in orbs:
                    orbs[-1].occupation = orbdata['occupations'][iorb]
            wfn.add_orbitals(orbs, orbtype=setname)

        return wfn

    @staticmethod
    def _json_to_quantities(jsonprops):  # TODO: handle this within JSON decoder
        props = {}
        for name, property in jsonprops.items():
            if isinstance(property, dict) and len(property) == 2 and \
                            'units' in property and 'value' in property:
                props[name] = property['value'] * u.ureg(property['units'])
            else:
                props[name] = property

        return props

    def _get_inputfiles(self):
        """ Override this method to pass additional input files to the program
        """
        return {}

    def minimize(self, nsteps=None):
        job = self._make_minimization_job(nsteps)

        return mdt.compute.run_job(job, _return_result=True)

    def _make_minimization_job(self, nsteps):
        params, inputfiles = self._prep_calculation([self.DEFAULT_PROPERTIES])
        params['runType'] = 'minimization'
        if nsteps is not None:
            params['minimization_steps'] = 100
        inputfiles['params.json'] = mdt.utils.json_dumps(dict(params))

        job = pyccc.Job(image=mdt.compute.get_image_path(self.IMAGE),
                        command='%s && %s' % (self.RUNNER, self.PARSER),
                        inputs=inputfiles,
                        when_finished=self.finish_min,
                        name='%s/%s' % (self.MODELNAME, self.mol.name))
        return job

    def finish_min(self, job):
        # TODO: parse more data than just the final minimization state
        traj = mdt.Trajectory(self.mol)
        traj.new_frame()

        results = json.loads(job.get_output('results.json').read())
        new_state = self._json_to_quantities(results['states'][0])

        self.mol.positions = new_state['positions']
        self.mol.properties = self._process_results(results)
        traj.new_frame()
        return traj

