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
import pyccc
import pyccc.exceptions

from moldesign import compute, exceptions
from moldesign.interfaces.openmm import force_remote, MdtReporter, pint2simtk, OpenMMPickleMixin
from moldesign.utils import exports


from .base import IntegratorBase, LangevinBase


class OpenMMBaseIntegrator(IntegratorBase, OpenMMPickleMixin):
    _openmm_compatible = True

    def prep(self):
        if not self.mol.energy_model._openmm_compatible:
            raise exceptions.NotSupportedError(
                    ("%s requires an OpenMM-based energy model (current "
                     "energy model '%s' is not supported with this integrator") %
                    (self.mol.energy_model.__class__.__name__, self.__class__.__name__))

        # For OpenMM integrators, all system setup is delegated to the energy model
        # this will create values for self.energy_model, self.sim, and and self._prepped
        self.mol.energy_model.prep()

    def run(self, run_for, wait=False):
        try:
            traj = self._run(run_for)
        except pyccc.exceptions.ProgramFailure:
            raise pyccc.exceptions.ProgramFailure(
                    'OpenMM crashed silently. Please examine the output. '
                    'This may be due to large forces from, for example, '
                    'an insufficiently minimized starting geometry.')
        if force_remote or (not wait):
            self.mol.energy_model._sync_remote(traj.mol)
            traj.mol = self.mol
        return traj

    @compute.runsremotely(enable=force_remote, is_imethod=True)
    def _run(self, run_for):
        self.prep()
        self.energy_model._set_constraints()  # calling this just to raise an exception if necessary
        nsteps = self.time_to_steps(run_for, self.params.timestep)
        self.energy_model._set_openmm_state()

        self.reporter = self._attach_reporters()
        self.reporter.annotation = self._describe_dynamics()
        self.reporter.report_from_mol()

        self.sim.step(nsteps)  # this is the actual dynamics loop

        self.energy_model._sync_to_openmm()
        if self.reporter.last_report_time != self.mol.time:
            self.reporter.report_from_mol()
        return self.reporter.trajectory

    def _attach_reporters(self):
        """
        Make sure the simulation has reporters for this run
        :return:
        """
        report_interval = self.time_to_steps(self.params.frame_interval,
                                             self.params.timestep)
        reporter = MdtReporter(self.mol, report_interval)
        self.sim.reporters = [reporter]
        return reporter

    def _describe_dynamics(self):
        raise NotImplementedError()


@exports
class OpenMMVerlet(OpenMMBaseIntegrator):
    def get_openmm_integrator(self):
        from simtk import openmm
        integrator = openmm.VerletIntegrator(pint2simtk(self.params.timestep))
        return integrator

    def _describe_dynamics(self):
        return 'Constant energy dynamics'


@exports
class OpenMMLangevin(LangevinBase, OpenMMBaseIntegrator):
    def get_openmm_integrator(self):
        from simtk import openmm

        integrator = openmm.LangevinIntegrator(
            pint2simtk(self.params.temperature),
            pint2simtk(self.params.collision_rate),
            pint2simtk(self.params.timestep))
        return integrator

    def _describe_dynamics(self):
        return 'Langevin dynamics @ %s' % self.params.temperature



