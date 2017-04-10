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
class LAMMPSLangevin(ConstantTemperatureBase):

    NAME_RESULT = "result"  # Name used for dump
    NAME_AFFECTED_ATOMS = "affected_atoms"  # Name used for grouped atoms
    NAME_LANGEVIN = "langevin_sim"    # Name used for Langevin
    NAME_ADDFORCE = "add_user_force"    # Name used for addforce
    
    # TODO: raise exception if any constraints are requested ...
    def __init__(self, *args, **kwargs):
        super(LAMMPSLangevin, self).__init__(*args, **kwargs)
        self._prepped = False # is the model prepped?
        self._model = None

        self.lammps_system = None
        self.total_iter = None      # How many iterations of timesteps before end of simulation
        self.output_iter = None     # How often we dump the positions of the atoms during the simulation
        self.traj = None

    def run(self, run_for):
        """
        Users won't call this directly - instead, use mol.run
        Propagate position, momentum by a single timestep using LAMMPS Langevin
        
        Arguments:
            run_for (int): number of timesteps OR amount of time to run for
            self.params.timestep (float): timestep length

        """
        self.total_iter = self.time_to_steps(run_for, self.params.timestep)
        self.output_iter = self.time_to_steps(self.params.frame_interval, self.params.timestep)

        if self.total_iter < self.output_iter:
            raise ValueError("Run duration {0} "
                             "can\'t be smaller than frame interval {1}".format(self.total_iter, self.output_iter))

        # Prepare lammps system
        self.prep()
        self._configure_system()
        lmps = self.lammps_system

        predictable_filename = "dump.result"
        dump_filepath = os.path.join(self.mol.energy_model.tmpdir, predictable_filename)
        self.mol.energy_model.files.append(dump_filepath)

        dump_command = "dump {0} all custom {1} {2} id x y z vx vy vz".format(self.NAME_RESULT,
                                                                            self.output_iter, dump_filepath)

        lmps.command(dump_command)
        lmps.command("dump_modify {0} sort id".format(self.NAME_RESULT))

        # run simulation
        lmps.run(self.total_iter, "post no")

        # Reset lammps system by undumping and unfixing
        for dmp in lmps.dumps:
            lmps.command("undump " + dmp.get('name'))

        for fix in lmps.fixes:
            lmps.command("unfix " + fix.get('name'))

        # Create trajectory object
        self.traj = Trajectory(self.mol)
        self.mol.calculate()

        # Dynamics loop over the dump file
        dump = open(dump_filepath, "rw+")
        content = dump.readlines()
        for istep in xrange(self.total_iter/self.output_iter):
            self.step(istep, content)
            self.traj.new_frame()
        dump.close()

        # Set time
        self.mol.time += self.time

        return self.traj

    def prep(self):
        """
        Prepare the LAMMPS system by configuring it for Langevin simulation
        
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
        # Get lammps object from model
        lmps = self._model.lammps_system

        # Check if we need to constrain hbonds
        if self.params.get('constrain_hbonds', False):
            raise NotImplementedError('SHAKE not implemented')

        # Check if we need to constrain water
        if self.params.get('constrain_water', False):
            raise NotImplementedError('SHAKE not implemented')

        # Set timestep
        lmps.command("timestep " + str(self.params.timestep.value_in(u.fs)))
        
        # Set Langevin settings
        lmps.command("fix lange_nve all nve")
        langevin_command = "fix {0} all langevin {1} {2} {3} 48279".format(self.NAME_LANGEVIN,
                                                                self.params.temperature.value_in(u.kelvin),
                                                                self.params.temperature.value_in(u.kelvin), 100.0)
        lmps.command(langevin_command)

        # Group affected atom
        if self._model.affected_atoms is not None:
            group_atom_cmd = "group {0} id ".format(self.NAME_AFFECTED_ATOMS)
            atoms = self._model.affected_atoms
            for atom in atoms:
                group_atom_cmd += "{0} ".format(atom.index+1)
            print group_atom_cmd
            lmps.command(group_atom_cmd.rstrip())

        # Apply user force fix
        if self._model.force_vector is not None:
            force_vector = self._model.force_vector
            lmps.command("fix {0} all addforce {1} {2} {3}".format(self.NAME_ADDFORCE, force_vector[0],
                                                                           force_vector[1], force_vector[2]))
        lmps.command("thermo_style custom step")
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