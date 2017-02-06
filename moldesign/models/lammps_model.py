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
    """Creates an OpenMM "context" to drive energy calculations.
    Note that, while a dummy integrator is assigned, a different context will
    be created for any MD calculations.

    :ivar sim: openmm simulation object
    :type sim: simtk.openmm.app.Simulation
    """
    # NEWFEATURE: need to set/get platform (and properties, e.g. number of threads)
    DEFAULT_PROPERTIES = ['potential_energy', 'forces']
    
    PARAMETERS = [
        mdt.parameters.Parameter(
            'skin', 'Extra distance beyond force cutoff',
            type=u.angstrom, default=2.0)]

    
    def __init__(self, **kwargs):
        super(LAMMPSPotential, self).__init__(**kwargs)
        self._prepped = False # is the model prepped?

        self.lammps_system = None # LAMMPS object for this molecule
        self.hydrogen_atom_types = None # are water molecules grouped?
        self._last_velocity = None # keep track of previous velocities of the molecule
        self.unit_system = None


    # calculate potential energy and force
    def calculate(self, requests):
        
        # Recreate system if molecule changed
        self.prep()

        # Run for 0 fs to evaluate system
        my_lmps = self.lammps_system
        my_lmps.run(0)

        # Evaluate potential energy and forces
        pe = my_lmps.eval("pe")
        force_array = []
        for i in range(0, my_lmps.atoms.natoms):
            force_array.append(my_lmps.atoms[i].force)
        return {'potential_energy': pe * u.kcalpermol,
                'forces': force_array * u.kcalpermol / u.angstrom}
    
      
    def prep(self):
        """
        Drive the construction of the LAMMPS simulation
        This will rebuild this OpenMM simulation if: A) it's not built yet, or B)
        there's a new integrator
        """ 
        
        # If current molecule's velocity is not the same as last recorded velocity, create a new system
        if self._last_velocity is None or (self._last_velocity == self.mol.velocities).all():
            self._create_system()

        self.unit_system = u.UnitSystem(length=u.angstrom, force=u.kcalpermol/u.fs, energy=u.kcalpermol, time=u.fs, 
            mass=u.gpermol, charge=u.electron_charge)

        self._prepped = True
        
        # TODO: wiring integrator

    
    #################################################
    # "Private" methods for managing LAMMPS are below
    
    def _create_system(self):
        """
        Create a LAMMPS system. Use MDT molecule object to construct a LAMMPS system
        
        Arguments:
            run_for (int): number of timesteps OR amount of time to run for
            self.params.timestep (float): timestep length
        """   
        # Ensure force fields are assigned 
        if 'amber_params' not in self.mol.ff:
            raise NotImplementedError('Assign forcefield parameters to the system')

        # Store the PRMTOP file if ff are assigned and use it to create lammps molecule data
        self.mol.ff.amber_params.prmtop.put('/tmp/tmp.prmtop')

        parmedmol = med.load_file('/tmp/tmp.prmtop')

        # initialize lammps object
        lmp = lammps()
        pylmp = PyLammps(ptr=lmp)

        # create temporary file system
        tmpdir = tempfile.mkdtemp()
        saved_umask = os.umask(0077)
        dataPath = self.create_lammps_data(tmpdir, parmedmol)

        pylmp.command("units real")
        pylmp.command("dimension 3")
        pylmp.command("atom_style full")
        pylmp.command("pair_style lj/charmm/coul/long 8.0 10.0 10.0")
        pylmp.command("bond_style harmonic")
        pylmp.command("angle_style harmonic")
        pylmp.command("dihedral_style harmonic")
        pylmp.command("improper_style harmonic")
        pylmp.command("kspace_style pppm 0.0001")
        pylmp.command("read_data " + dataPath)

        pylmp.command("neighbor    1 multi")
        pylmp.command("neigh_modify    delay 0")

        # securely remove temporary filesystem
        os.remove(dataPath)
        os.umask(saved_umask)
        os.rmdir(tmpdir)
        
        # group hbonds
        self.group_hydrogen(parmedmol)

        self.lammps_system = pylmp
        
        self._last_velocity = self.mol.velocities.copy() #Keep track of last velocity that was used to create the LAMMPS system 

    
    def group_hydrogen(self, parmedmol):
        """
        Get indices of non-bond atom types that contain hydrogen bonds
        
        Arguments:
            parmed: parmed struct to iterate through all nonbond types to find hbonds

        """
        if len(parmedmol.LJ_types) > 0:    
            for nonbond_name, nonbond_idx in parmedmol.LJ_types.iteritems():
                if(nonbond_name == "H"):
                    self.hydrogen_atom_types = nonbond_idx
                    break


    # In order to simluate using LAMMPS, we need:
    #   get molecule from pdb, cif, or name
    #   assign forcefields
    #   load parmed structure


    # LAMMPS UNITS:
    # For style real, these are the units:

    # mass = grams/mole
    # distance = Angstroms
    # time = femtoseconds
    # energy = Kcal/mole
    # velocity = Angstroms/femtosecond
    # force = Kcal/mole-Angstrom
    # torque = Kcal/mole
    # temperature = Kelvin
    # pressure = atmospheres
    # dynamic viscosity = Poise
    # charge = multiple of electron charge (1.0 is a proton)
    # dipole = charge*Angstroms
    # electric field = volts/Angstrom
    # density = gram/cm^dim

    # See http://parmed.github.io/ParmEd/html/index.html for ParmEd units

    
    def create_lammps_data(self, tmpdir, parmedmol):
        """
        Create a data. file that contains all necessary information regarding the molecule
        in order to create the appropirate LAMMPS system for it
        
        Arguments:
            tmpdir: temporary directory to store the data file
            parmed: parmed struct to iterate through all nonbond types to find hbonds

        """    
        predictable_filename = 'data.lammps_mol'

        path = os.path.join(tmpdir, predictable_filename)

        datalines =  ""
        datalines += "LAMMPS Description\r\n\r\n"
        datalines += "{0} atoms\r\n".format(len(parmedmol.atoms))
        datalines += "{0} bonds\r\n".format(len(parmedmol.bonds))
        datalines += "{0} angles\r\n".format(len(parmedmol.angles))
        datalines += "{0} dihedrals\r\n".format(len(parmedmol.dihedrals))
        datalines += "{0} impropers\r\n".format(len(parmedmol.impropers))
        datalines += "\r\n"
                    
        # Get non-bond atom types
        parmedmol.fill_LJ()
        max_idx = 0
        for nonbond_idx in parmedmol.LJ_types.itervalues():
            if(nonbond_idx > max_idx):
                max_idx = nonbond_idx

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


        # calculate lo and hi coordinates for Box size
        # xlo = xhi = ylo = yhi = zlo = zhi = None
        # for atom in self.mol.atoms:
        #     if xlo == None:
        #         xlo = atom.x.value_in(u.angstrom)
        #         xhi = atom.x.value_in(u.angstrom)
        #         ylo = atom.y.value_in(u.angstrom)
        #         yhi = atom.y.value_in(u.angstrom)
        #         zlo = atom.z.value_in(u.angstrom)
        #         zhi = atom.z.value_in(u.angstrom)
        #     else:
        #         xlo = min(xlo, atom.x.value_in(u.angstrom))
        #         xhi = max(xhi, atom.x.value_in(u.angstrom))
        #         ylo = min(ylo, atom.y.value_in(u.angstrom))
        #         yhi = max(yhi, atom.y.value_in(u.angstrom))
        #         zlo = min(zlo, atom.z.value_in(u.angstrom))
        #         zhi = max(zhi, atom.z.value_in(u.angstrom))

        datalines += "{0} {1} xlo xhi\r\n".format(-50, 100)
        datalines += "{0} {1} ylo yhi\r\n".format(-50, 100)
        datalines += "{0} {1} zlo zhi\r\n".format(-50, 100)

        print datalines

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

        with open(path, "w") as lammps_data:   
            lammps_data.write(datalines)
            
       # print datalines

        return path


