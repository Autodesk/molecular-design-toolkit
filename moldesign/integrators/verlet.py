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
import moldesign as mdt
from moldesign import units as u
from moldesign.molecules import Trajectory

from .base import IntegratorBase


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class VelocityVerlet(IntegratorBase):
    def __init__(self, *args, **kwargs):
        super(VelocityVerlet, self).__init__(*args, **kwargs)

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
        for istep in xrange(nsteps):
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
