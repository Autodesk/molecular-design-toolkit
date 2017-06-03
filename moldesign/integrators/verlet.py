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
from .. import units as u
from ..molecules import Trajectory
from ..utils import exports

from .base import IntegratorBase



@exports
class VelocityVerlet(IntegratorBase):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    # TODO: raise exception if any constraints are requested ...

    def run(self, run_for):
        """
        Users won't call this directly - instead, use mol.run
        Propagate position, momentum by a single timestep using velocity verlet
        :param run_for: number of timesteps OR amount of time to run for
        """
        if not self._prepped:
            self.prep()
        nsteps = self.time_to_steps(run_for, self.params.timestep)

        # Set up trajectory and record the first frame
        self.mol.time = 0.0 * u.default.time
        self.traj = Trajectory(self.mol)
        self.mol.calculate()
        self.traj.new_frame()
        next_trajectory_frame = self.params.frame_interval

        # Dynamics loop
        for istep in range(nsteps):
            self.step()
            if istep + 1 >= next_trajectory_frame:
                self.traj.new_frame()
                next_trajectory_frame += self.params.frame_interval
        return self.traj

    def prep(self):
        self.time = 0.0 * self.params.timestep
        self._prepped = True

    def step(self):
        # Move momenta from t-dt to t-dt/2
        phalf = self.mol.momenta + 0.5 * self.params.timestep * self.mol.calc_forces(wait=True)

        # Move positions from t-dt to t
        self.mol.positions += phalf * self.params.timestep / self.mol.dim_masses

        # Move momenta from t-dt/2 to t - triggers recomputed forces
        self.mol.momenta = phalf + 0.5 * self.params.timestep * self.mol.calc_forces(wait=True)
        self.time += self.params.timestep
        self.mol.time = self.time
