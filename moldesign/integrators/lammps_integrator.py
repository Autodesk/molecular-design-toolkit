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
import numpy as np

from .base import IntegratorBase, ConstantTemperatureBase

try:
    from lammps import lammps, PyLammps
    # WARNING: this is the real library, not our interface - this works because of absolute
    # imports. We should probably rename this interface
except ImportError:
    print 'LAMMPS could not be imported; using remote docker container'
    force_remote = True
else:  # this should be configurable
    force_remote = False  # debugging

import tempfile
import os
import sys
from itertools import islice
import time

def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class LAMMPSNvt(ConstantTemperatureBase):
    
    # TODO: raise exception if any constraints are requested ...
    def __init__(self, *args, **kwargs):
        super(LAMMPSNvt, self).__init__(*args, **kwargs)
        self._prepped = False # is the model prepped?
        self.lammps_system = None        
        
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

        tot_it = self.time_to_steps(run_for, self.params.timestep)
        output_freq = self.time_to_steps(self.params.frame_interval, self.params.timestep)
        
        if tot_it < output_freq:
            raise ValueError("Run duration can\'t be smaller than frame interval")

        #int(self.params.frame_interval.value_in(u.fs) / self.params.timestep.value_in(u.fs))
        
        # create temporary file system
        tmpdir = tempfile.mkdtemp()
        saved_umask = os.umask(0077)

        pos_path = os.path.join(tmpdir, "dump.position")
        vel_path = os.path.join(tmpdir, "dump.velocity")

        lmps = self.lammps_system

        lmps.command("dump position all custom {0} {1} id x y z".format(output_freq, pos_path))
        lmps.command("dump velocity all custom {0} {1} id vx vy vz".format(output_freq, vel_path))
        lmps.command("dump_modify position sort id")
        lmps.command("dump_modify velocity sort id")
        lmps.command("thermo_style custom step temp pe etotal")
        lmps.command("thermo 5")
        # Set up trajectory and record the first frame
        self.mol.time = 0.0 * u.fs
        lammps_units = self.model.unit_system
        self.traj = Trajectory(self.mol, unit_system=lammps_units)
        self.mol.calculate()

        # run simulation
        lmps.run(tot_it)

        # Dynamics loop
        for istep in xrange(tot_it/output_freq):
            self.step(istep, pos_path, vel_path)
            self.traj.new_frame()

        lmps.command("undump position")
        lmps.command("undump velocity")
        self.lammps_system = lmps
         
        # securely remove temporary filesystem
        os.remove(pos_path)
        os.remove(vel_path)
        os.umask(saved_umask)
        os.rmdir(tmpdir)
        
        return self.traj

        

    
    def prep(self):
        """
        Prepare the LAMMPS system by configuring it for NVT simulation
        
        Arguments:
            run_for (int): number of timesteps OR amount of time to run for
            self.params.timestep (float): timestep length

        """
        # sanity check of parameters
        # NOTE: DO I NEED THIS??
        if self.params.timestep.value_in(u.fs) > self.params.frame_interval.value_in(u.fs):
            raise ValueError

        #if self._prepped and self.model is self.mol.energy_model and self.model._prepped: return
        
        # prepare lammps model prior to preparing the integrator
        self.model = self.mol.energy_model
        self.model.prep()

        self.time = 0.0 * self.params.timestep
        self._prepped = True

        # get lammps object from model
        lammps_system = self.model.lammps_system

        # NOTE: Ensure time step is in femtoseconds 
        print self.params.timestep.value_in(u.fs)
        lammps_system.command("timestep " + str(self.params.timestep.value_in(u.fs)))
        
        # TODO:
        nvt_command = "fix 1 all nvt temp {0} {1} {2} tchain 1" .format(self.params.temperature.value_in(u.kelvin), 
            self.params.temperature.value_in(u.kelvin), 100.0)
        lammps_system.command(nvt_command)

        # NOTE: Can you help me figure out the parameters for fix shake?
        # if self.params.constrain_hbonds and self.model.hbond_exists == True :
        #     shake_hbond_command = "fix 2 hbond rattle 0.0001 10 0 b 4 6 8 10 12 14 18 a 31"
        #     lammps_system.command(shake_hbond_command)

        # NOTE: Can you help me figure out the parameters for fix shake?
        # if self.params.constrain_water and self.model.water_exists == True :
        #     shake_water_command = "fix 3 water rattle 0.0001 10 0 b 4 6 8 10 12 14 18 a 31"
        #     lammps_system.command(shake_water_command)

        self.lammps_system = lammps_system

    """
        Update molecule's position and velocity at each step of the simulation

    """
    def step(self, idx, pos_filepath, vel_filepath):
        lmps = self.lammps_system
        
        po = open(pos_filepath, "rw+")
        # start at 9, stop at 9+L.atoms.natoms
        # start at 14+L.atoms.natoms
        tmp = list(islice(po.readlines(), 9*(idx+1) + lmps.atoms.natoms*idx, (lmps.atoms.natoms+9)*(idx+1)))
        tmp = list(item.split()[1:] for item in tmp)
        pos_array = np.array(tmp).astype(np.float)

        self.mol.positions = pos_array * u.angstrom
        # print self.mol.positions[1004] # Debugging
        # print pos_array
        # print "\r\n"
        vel = open(vel_filepath, "rw+")
        tmp = list(islice(vel.readlines(), 9*(idx+1) + lmps.atoms.natoms*idx, (lmps.atoms.natoms+9)*(idx+1)))
        tmp = list(item.split()[1:] for item in tmp)
        vel_array = np.array(tmp).astype(np.float)

        self.mol.velocities = vel_array * u.angstrom / u.fs
        # print self.mol.velocities[1004] # Debugging
        self.time += self.params.frame_interval
        self.mol.time = self.time

    #################################################
    # "Private" methods for managing LAMMPS are below





