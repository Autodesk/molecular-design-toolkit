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
import moldesign.units as u

from moldesign import compute
from moldesign import forcefields as ff

from .base import EnergyModelBase

try:
    from lammps import lammps, PyLammps
    # WARNING: this is the real library, not our interface - this works because of absolute
    # imports. We should probably rename this interface
except ImportError:
    print 'LAMMPS could not be imported; using remote docker container'
    force_remote = True
else:  # this should be configurable
    force_remote = False  # debugging

debugging = True # debugging

import parmed as med
import tempfile
import os
import sys
from .amber2lammps import *

def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class LAMMPSPotential(EnergyModelBase):
    """Creates a LAMMPS system to drive energy calculations.
    """
    # NEWFEATURE: need to set/get platform (and properties, e.g. number of threads)
    DEFAULT_PROPERTIES = ['potential_energy', 'forces']

    def __init__(self, **kwargs):
        super(LAMMPSPotential, self).__init__(**kwargs)
        self._prepped = False # is the model prepped?
        self._lammps_data_path = None

        # For temporary directory
        # create temporary file system
        # Create temporary filesystem
        self.tmpdir = tempfile.mkdtemp()
        self.saved_umask = os.umask(0077)
        self.files = []

        self.lammps_system = None # Basic LAMMPS system for this molecule
        self.unit_system = None

    def __exit__(self, exc_type, exc_value, traceback):
        # Unlink all files
        for f in self.files:
            os.unlink(f)
            os.remove(f)

        # securely remove temporary file system
        os.umask(self.saved_umask)
        os.rmdir(self.tmpdir)

        # Delete lammps object
        del self.lammps_system

    def prep(self):
        """
        Drive the construction of the LAMMPS simulation
        This will rebuild this OpenMM simulation if: A) it's not built yet, or B)
        there's a new integrator
        """

        if not self._prepped or self.lammps_system is None:
            # Create lammps
            self._create_lammps()
            self.unit_system = u.UnitSystem(length=u.angstrom, force=u.kcalpermol / u.fs, energy=u.kcalpermol,
                                            time=u.fs,
                                            mass=u.gpermol, charge=u.electron_charge)

        self._prepped = True

    # calculate potential energy and force
    def calculate(self, requests):
        
        # Recreate system if molecule changed
        results = self._evaluate_molecule()
        return results

    #################################################
    # "Private" methods for managing LAMMPS are below

    def _evaluate_molecule(self):
        # Recreate system if molecule changed
        self.prep()

        # Run for 0 fs to evaluate system
        my_lmps = self.lammps_system
        my_lmps.command("thermo_style one")
        my_lmps.run(0)

        # Evaluate potential energy and forces
        pe = my_lmps.eval("pe")
        force_array = []
        for i in range(0, my_lmps.atoms.natoms):
            force_array.append(my_lmps.atoms[i].force)
        return {'potential_energy': pe * u.kcalpermol,
                'forces': force_array * u.kcalpermol / u.angstrom}

    # NOTE: Used by workflows to create lammps formatted data file
    def get_lammps_data_path(self):
        self.prep()

        # initialize lammps object
        return self._lammps_data_path
    
    def _create_lammps(self):

        """
        Create a LAMMPS system. Use MDT molecule object to construct a LAMMPS system

        """
        crd_filepath = '/tmp/tmp.crd'
        top_filepath = '/tmp/tmp.top'

        # Put amber parameters
        self.mol.ff.amber_params.prmtop.put(top_filepath)
        self.mol.ff.amber_params.inpcrd.put(crd_filepath)

        # Add file paths to list of files
        self.files.append(crd_filepath)
        self.files.append(top_filepath)

        # Ensure force fields are assigned 
        if 'amber_params' not in self.mol.ff:
            raise NotImplementedError('Assign forcefield parameters to the system')

        # Store the PRMTOP file if ff are assigned and use it to create lammps molecule data

        # initialize lammps object
        lmp = lammps()
        pylmp = PyLammps(ptr=lmp)

        # Create temperary data file
        predictable_filename = 'data.lammps_mol'
        data_path = os.path.join(self.tmpdir, predictable_filename)

        # Add file paths to list of files
        self.files.append(data_path)
        self._lammps_data_path = data_path

        # Convert amber params files to lammps data
        Convert_Amber_files(crd_filepath, top_filepath, data_path)

        # Set basic parameters. Convert_Amber_files uses harmonic style
        pylmp.command("units real")
        pylmp.command("dimension 3")
        pylmp.command("atom_style full")
        pylmp.command("pair_style lj/charmm/coul/long 8.0 10.0 10.0")
        pylmp.command("bond_style harmonic")
        pylmp.command("angle_style harmonic")
        pylmp.command("dihedral_style harmonic")
        pylmp.command("improper_style harmonic")
        pylmp.command("kspace_style pppm 0.0001")
        pylmp.command("read_data " + data_path)

        # Assign lammps system
        self.lammps_system = pylmp

@exports
class LAMMPSInteractive(EnergyModelBase):

    def __init__(self, **kwargs):
        super(LAMMPSInteractive, self).__init__(**kwargs)
        self._prepped = False   # is the model prepped?

        # LAMMPSPotential model object
        self._potential_model = None
        self.lammps_system = None  # Get from potential model
        self.unit_system = None    # Get from potential model
        self.tmpdir = None
        self.saved_umask = None
        self.files = None

        self._calculated_properties = None
        self.affected_atoms = None  # atom that the user is directly interacting with
        self.force_vector = None   # force that acts on the molecule

    def __exit__(self, exc_type, exc_value, traceback):
        # Delete LAMMPS potential model
        del self._potential_model

    def prep(self):
        """
        Drive the construction of the LAMMPS simulation
        This will rebuild this OpenMM simulation if: A) it's not built yet, or B)
        there's a new integrator
        """

        if not self._prepped or self._potential_model is None or self.lammps_system is None:
            self._potential_model = mdt.models.LAMMPSPotential()
            self._potential_model.mol = self.mol

            self._potential_model.prep()

            self._calculated_properties = None

            # Get settings from LAMMPSPotential model
            self.lammps_system = self._potential_model.lammps_system
            self.unit_system = self._potential_model.unit_system
            self.tmpdir = self._potential_model.tmpdir
            self.saved_umask = self._potential_model.saved_umask
            self.files = self._potential_model.files

        self._prepped = True

    def calculate(self, requests):
        '''
            After assigning this energy model to the molecule, call this function as soon as possible, so
            the molecule is ready for the interactive simulation
        '''
        
        # Recreate system if molecule changed
        if self._calculated_properties is None:
            self._calculated_properties = self._potential_model.calculate()

        return self._calculated_properties
    
    def apply_force(self, atoms, vector):
        '''
            :param atoms: atoms to apply force
            :param vector: amount of force to apply to atoms in kcal/mol/angstrom
        '''
        if atoms is None or vector is None:
            raise ValueError('One or more argument(s) missing!')

        vector_real_units = vector.value_in(u.kcalpermol / u.angstrom)
        
        # Set internal variables
        self.affected_atoms = atoms
        self.force_vector = vector_real_units

    def reset_force(self):
        '''
            Reset user-added force
        '''
        self.affected_atoms = None
        self.force_vector = None

    def get_lammps_data(self):
        self.prep()
        lammps_data_path = self._potential_model.get_lammps_data_path()

        export_data_string = ''
        with open(lammps_data_path) as data:
            lines = data.readlines()
            for line in lines:
                export_data_string += line
        data.close()

        return export_data_string



