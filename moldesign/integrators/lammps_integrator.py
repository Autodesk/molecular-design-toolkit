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
from __future__ import absolute_import

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
from copy import deepcopy
from itertools import islice

def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class LAMMPSNvt(ConstantTemperatureBase):

    NAME_RESULT = "result"
    NAME_NVT_SIM = "nvt_sim"
    
    # TODO: raise exception if any constraints are requested ...
    def __init__(self, *args, **kwargs):
        super(LAMMPSNvt, self).__init__(*args, **kwargs)
        self._prepped = False # is the model prepped?
        self._model = None
        self.lammps_system = None
        self.traj = None
        
    def run(self, run_for):
        """
        Users won't call this directly - instead, use mol.run
        Propagate position, momentum by a single timestep using LAMMPS NVT
        
        Arguments:
            run_for (int): number of timesteps OR amount of time to run for
            self.params.timestep (float): timestep length

        """

        # Prepare lammps system
        self.prep()
        self._configure_system()

        tot_it = self.time_to_steps(run_for, self.params.timestep)
        if type(self._model) is mdt.models.LAMMPSInteractive:
            self.params.frame_interval = run_for
            output_freq = tot_it
        else:
            output_freq = self.time_to_steps(self.params.frame_interval, self.params.timestep)
        
        if tot_it < output_freq:
            raise ValueError("Run duration {0} can\'t be smaller than frame interval {1}".format(tot_it, output_freq))
        
        # create temporary file system
        tmpdir = tempfile.mkdtemp()
        saved_umask = os.umask(0077)
        result_path = os.path.join(tmpdir, "dump.result")

        lmps = self.lammps_system

        # TODO: shake not working properly
        # if self.params.constrain_hbonds and self._model.hbond_group is not None:
        #     shake_cmd = "fix constrain_hbonds all shake 0.0001 50 0 m {0}".format(self._model.hbond_group)
        #     # print shake_cmd
        #     lmps.command(shake_cmd)
        #     lmps.command("run_style respa 2 2")

        lmps.command("dump {0} all custom {1} {2} id x y z vx vy vz".format(self.NAME_RESULT, output_freq, result_path))
        lmps.command("dump_modify {0} sort id".format(self.NAME_RESULT))

        # run simulation
        lmps.run(tot_it, "post no")

        # Reset lammps system by undumping and unfixing
        for dmp in lmps.dumps:
            lmps.command("undump " + dmp.get('name'))

        for fix in lmps.fixes:
            lmps.command("unfix " + fix.get('name'))

        # Read dynamic simulation results
        dump = open(result_path, "rw+")
        results = dump.readlines()

        # Create trajectory unless interactive model
        if type(self._model) is not mdt.models.LAMMPSInteractive:
            self.traj = Trajectory(self.mol, unit_system=self._model.unit_system)
            self.mol.calculate()

        # Dynamics loop
        for istep in xrange(tot_it/output_freq):
            # Start reading from index 1 since the first round of positions/vectors 
            # would be the same as the positions before the simulation
            self.step(istep+1, results)
            if type(self._model) is not mdt.models.LAMMPSInteractive:
                self.traj.new_frame()

        self.mol.time += self.time

        # securely remove temporary filesystem
        os.remove(result_path)
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
        # prepare lammps model prior to preparing the integrator
        self._model = self.mol.energy_model
        self._model.prep()

        self.time = 0.0 * self.params.timestep
        self._prepped = True

    def _configure_system(self):
        # get lammps object from model
        lmps = self._model.lammps_system

        # NOTE: Ensure time step is in femtoseconds
        lmps.command("timestep " + str(self.params.timestep.value_in(u.fs)))
        # print self.params.timestep.value_in(u.fs)

        nvt_command = "fix {0} all nvt temp {1} {2} {3}".format(self.NAME_NVT_SIM,
                                                                self.params.temperature.value_in(u.kelvin),
                                                                self.params.temperature.value_in(u.kelvin), 100.0)
        # print nvt_command
        lmps.command(nvt_command)

        lmps.command("thermo_style custom step")
        # lmps.command("thermo 1000")

        self.lammps_system = lmps  # Save lammps configuration

    def step(self, idx, results):
        
        """
            Update molecule's position and velocity at each step of the simulation

        """
        lmps = self.lammps_system
        # start at 9, stop at 9+L.atoms.natoms
        # start at 14+L.atoms.natoms
        
        tmp = list(islice(results, 9*(idx+1) + lmps.atoms.natoms*idx, (lmps.atoms.natoms+9)*(idx+1)))
        pos = list(item.split()[1:4] for item in tmp)
        vel = list(item.split()[4:] for item in tmp)

        pos_array = np.array(pos).astype(np.float)
        vel_array = np.array(vel).astype(np.float)
        
        self.mol.positions = pos_array * u.angstrom
        self.mol.velocities = vel_array * u.angstrom / u.fs

        self.time += self.params.frame_interval





