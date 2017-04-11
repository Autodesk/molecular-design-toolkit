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

    # TODO:
    # def __del__(self, exc_type, exc_value, traceback):
    #     print "Deleting LAMMPSPotential"

    #     # Unlink all files
    #     for f in self.files:
    #         os.unlink(f)
    #         os.remove(f)

    #     # securely remove temporary file system
    #     os.umask(self.saved_umask)
    #     os.rmdir(self.tmpdir)

    #     # Delete lammps object
    #     del self.lammps_system

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
    def calculate(self, requests=None):
        
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
    
        # Ensure force fields are assigned 
        if self.mol.ff is None or self.mol.ff.parmed_obj is None:
            raise NotImplementedError('Assign forcefield parameters to the system')

        # initialize lammps object
        lmp = lammps()
        pylmp = PyLammps(ptr=lmp)

        # Create temperary data file
        predictable_filename = 'data.lammps_mol'
        data_path = os.path.join(self.tmpdir, predictable_filename)

        # Add file paths to list of files
        self.files.append(data_path)
        self._lammps_data_path = data_path

        with open(data_path, "w") as lammps_data:
            formatted_data = self._get_lammps_data()
            lammps_data.write(formatted_data)

        # Set basic parameters.
        pylmp.command("units real")
        pylmp.command("dimension 3")
        pylmp.command("atom_style full")
        pylmp.command("pair_style lj/charmm/coul/charmm/implicit 8.0 10.0")

        pylmp.command("bond_style harmonic")
        pylmp.command("angle_style harmonic")
        pylmp.command("dihedral_style harmonic")
        pylmp.command("improper_style harmonic")
        
        pylmp.command("dielectric 4.0")
        pylmp.command("read_data " + data_path)

        # Assign lammps system
        self.lammps_system = pylmp

    def _get_lammps_data(self):
        """
        Parse molecule and force fields info to generate a formatted data

        """    

        parmedmol = self.mol.ff.parmed_obj

         # Get non-bond atom types
        parmedmol.fill_LJ()
        max_idx = 0
        for nonbond_idx in parmedmol.LJ_types.itervalues():
            if nonbond_idx > max_idx:
                max_idx = nonbond_idx

        datalines =  ""
        
        datalines += "LAMMPS Description\r\n\r\n"
        datalines += "{0} atoms\r\n".format(len(parmedmol.atoms))
        datalines += "{0} bonds\r\n".format(len(parmedmol.bonds))
        datalines += "{0} angles\r\n".format(len(parmedmol.angles))
        datalines += "{0} dihedrals\r\n".format(len(parmedmol.dihedrals))
        datalines += "{0} impropers\r\n".format(len(parmedmol.impropers))
        datalines += "\r\n"
                         
        datalines += "{0} atom types\r\n".format(max_idx)
        if len(parmedmol.bond_types) > 0:
            datalines += "{0} bond types\r\n".format(len(parmedmol.bond_types))
        if len(parmedmol.angle_types) > 0:
            datalines += "{0} angle types\r\n".format(len(parmedmol.angle_types))
        if len(parmedmol.dihedral_types) > 0:
            datalines += "{0} dihedral types\r\n".format(len(parmedmol.dihedral_types))
        if len(parmedmol.improper_types) > 0:
            datalines += "{0} improper types\r\n".format(len(parmedmol.improper_types))
    
        datalines += "\r\n"

        # TODO: calculate lo and hi coordinates for Box size
        xlo = xhi = ylo = yhi = zlo = zhi = None
        for atom in self.mol.atoms:
            if xlo == None:
                xlo = atom.x.value_in(u.angstrom)
                xhi = atom.x.value_in(u.angstrom)
                ylo = atom.y.value_in(u.angstrom)
                yhi = atom.y.value_in(u.angstrom)
                zlo = atom.z.value_in(u.angstrom)
                zhi = atom.z.value_in(u.angstrom)
            else:
                xlo = min(xlo, atom.x.value_in(u.angstrom))
                xhi = max(xhi, atom.x.value_in(u.angstrom))
                ylo = min(ylo, atom.y.value_in(u.angstrom))
                yhi = max(yhi, atom.y.value_in(u.angstrom))
                zlo = min(zlo, atom.z.value_in(u.angstrom))
                zhi = max(zhi, atom.z.value_in(u.angstrom))

        datalines += "{0} {1} xlo xhi\r\n".format(xlo-50, xhi+50)
        datalines += "{0} {1} ylo yhi\r\n".format(ylo-50, yhi+50)
        datalines += "{0} {1} zlo zhi\r\n".format(zlo-50, zhi+50)

        datalines += "\r\n"

        # Masses - atom types
        datalines += "Masses\r\n\r\n"
        for i in range(1, max_idx+1):
            for atom in parmedmol.atoms:
                if atom.nb_idx == i:
                    datalines += "{0} {1}\r\n".format(i, atom.mass)
                    break
        datalines += "\r\n"
        
        # Pair coefficients                   
        datalines += "Pair Coeffs\r\n\r\n"
        for i in range(1, max_idx+1):
            for atom in parmedmol.atoms:
                if atom.nb_idx == i:
                    datalines += "{0} {1} {2} {3} {4}\r\n".format(i, atom.epsilon, atom.sigma, atom.epsilon_14, atom.sigma_14)
                    break
        datalines += "\r\n"
                    
        # Bond Coeffs - Harmonic

        if parmedmol.bond_types:
            datalines += "Bond Coeffs\r\n\r\n"
            for bt in parmedmol.bond_types:
                datalines += "{0} {1} {2}\r\n".format(bt.idx+1, bt.k, bt.req)
            
            datalines += "\r\n"

        # Angle Coeffs - Harmonic
        if parmedmol.angle_types:
            datalines += "Angle Coeffs\r\n\r\n"
            for at in parmedmol.angle_types:
                datalines +="{0} {1} {2}\r\n".format(at.idx+1, at.k, at.theteq)
           
            datalines +="\r\n"


        # Dihedral Coeffs - Charmm
        if parmedmol.dihedral_types:
            datalines += "Dihedral Coeffs\r\n\r\n"
            for dt in parmedmol.dihedral_types:
                datalines += "{0} {1} {2} {3}\r\n".format(dt.idx+1, dt.phi_k, (-1 if dt.phase >= 180 else 1), dt.per)
            
            datalines += "\r\n" 


        # Improper Coeffs - Harmonic
        if parmedmol.improper_types:
            datalines += "Improper Coeffs\r\n\r\n"
            for it in parmedmol.improper_types:
                datalines += "{0} {1} {2}\r\n".format(it.idx+1, it.psi_k, it.psi_eq)
            
            datalines += "\r\n"


        # print atoms
        datalines += "Atoms\r\n\r\n"
        for atom in parmedmol.atoms:
            datalines += "{0} {1} {2} {3} ".format(atom.idx+1, atom.residue.idx+1, atom.nb_idx, atom.charge)
            if not self.mol.positions.any():        
                datalines += "{0} {1} {2} ".format(0.00, 0.00, 0.00)    
            else:
                pos = self.mol.positions[atom.idx]
                datalines += "{0} {1} {2} ".format(pos[0].value_in(u.angstrom), pos[1].value_in(u.angstrom), 
                    pos[2].value_in(u.angstrom))
            datalines += "{0} {1} {2}\r\n".format(0, 0, 0) # nx, ny, nz not defined    

        datalines += "\r\n"


        # print velocities
        datalines += "Velocities\r\n\r\n"
        for atom in self.mol.atoms:
            xvel = atom.vx.value_in(u.angstrom/u.fs)
            yvel = atom.vy.value_in(u.angstrom/u.fs)
            zvel = atom.vz.value_in(u.angstrom/u.fs)
            datalines += "{0} {1} {2} {3}\r\n".format(atom.index+1, xvel, yvel, zvel)
            
        datalines += "\r\n"


        # print bonds
        datalines += "Bonds\r\n\r\n"
        for idx, bond in enumerate(parmedmol.bonds):
            datalines += "{0} {1} {2} {3}\r\n".format(idx+1, bond.type.idx+1, 
              bond.atom1.idx+1, bond.atom2.idx+1)
            
        datalines += "\r\n"


        # print angles
        datalines += "Angles\r\n\r\n"
        for idx, angle in enumerate(parmedmol.angles):
            datalines += "{0} {1} {2} {3} {4}\r\n".format(idx+1, angle.type.idx+1,
              angle.atom1.idx+1, 
              angle.atom2.idx+1, angle.atom3.idx+1)
            
        datalines += "\r\n"


        # print dihedrals
        if parmedmol.dihedrals:
            datalines += "Dihedrals\r\n\r\n"
            for idx, di in enumerate(parmedmol.dihedrals):
                datalines += "{0} {1} {2} {3} {4} {5}\r\n".format(idx+1, di.type.idx+1, 
                  di.atom1.idx+1, di.atom2.idx+1, 
                  di.atom3.idx+1, di.atom4.idx+1)
                
            datalines += "\r\n"


        # print impropers
        if parmedmol.impropers:
            datalines += "Impropers\r\n\r\n"
            for im in parmedmol.impropers:
                datalines += "{0} {1} {2} {3} {4} {5}\r\n".format(im.idx+1, im.type.idx+1, 
                  im.atom1.idx+1, im.atom2.idx+1, 
                  im.atom3.idx+1, im.atom4.idx+1)
            
            datalines += "\r\n"

        # Need both atom types and LAMMPS formatted data for web
        return datalines

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

    # TODO:
    # def __del__(self, exc_type, exc_value, traceback):
    #     print "Deleting LAMMPSInteractive"

    #     # Delete LAMMPS potential model
    #     del self._potential_model

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

        self.prep()

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



