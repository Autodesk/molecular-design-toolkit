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
class LAMMPSNvt(ConstantTemperatureBase):
    
    # TODO: raise exception if any constraints are requested ...

    def run(self, run_for):
        """
        Users won't call this directly - instead, use mol.run
        Propagate position, momentum by a single timestep using LAMMPS NVT
        
        Arguments:
            run_for (int): number of timesteps OR amount of time to run for
            self.params.timestep (float): timestep length

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

    """
        Prepare the LAMMPS system by configuring it for NVT simulation
        
        Arguments:
            run_for (int): number of timesteps OR amount of time to run for
            self.params.timestep (float): timestep length

    """
    def prep(self):

        # sanity check of parameters
        # NOTE: DO I NEED THIS??
        if self.params.timestep > params.frame_interval:
            raise ValueError

        if self._prepped and self.model is self.mol.energy_model and self.model._prepped: return
        
        # prepare lammps model prior to preparing the integrator
        self.model = self.mol.energy_model
        self.model.prep()

        self.time = 0.0 * self.params.timestep
        self._prepped = True

        # get lammps object from model
        lammps_system = self.model.lammps_system

        # NOTE: Ensure time step is in femtoseconds 
        lammps_system.command("timestep " + str(self.params.timestep))
        lammps_system.command("thermo_style custom step temp pe etotal")
        lammps_system.command("thermo " + str(self.params.frame_interval))
        
        # TODO:
        nvt_command = "fix 1 all nvt temp {0} {1} {2}" .format(self.params.temperature, 
            self.params.temperature, 100.0)
        lammps_system.command(nvt_command)

        # NOTE: Can you help me figure out the parameters for fix shake?
        if self.params.constrain_hbonds and self.model.group_hbond == True :
            shake_hbond_command = "fix 2 hbond shake 0.0001 20 0 t 5 6 m 1.0 a 31"
            lammps_system.command(shake_hbond_command)

        # NOTE: Can you help me figure out the parameters for fix shake?
        if self.params.constrain_water and self.model.group_water == True :
            shake_water_command = "fix 2 water shake 0.0001 20 0 t 5 6 m 1.0 a 31"
            lammps_system.command(shake_water_command)

        self.lammps_system = lammps_system


    """
        Update molecule's position and velocity at each step of the simulation

    """
    def step(self):
        # Run Lammps simulation
        L = self.lammps_system
        L.run(1)    # run takes in integer number

        # Update position and velocity of each atom
        for i in range(0, L.atoms.natom):
            self.mol.positions[i] = L.atoms[i].position
            self.mol.velocities[i] = L.atoms[i].velocity

        self.time += self.params.timestep
        self.mol.time = self.time

    #################################################
    # "Private" methods for managing LAMMPS are below





