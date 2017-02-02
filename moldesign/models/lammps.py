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
        self.group_hbond = False # are hbonds grouped?
        self.group_water = False # are water molecules grouped?
        self._last_velocity = None # keep track of previous velocities of the molecule
        print "Initialized!"


    # calculate potential energy and force
    def calculate(self, requests=None):
        
        # Recreate system if molecule changed
        # NOTE: WILL THIS WORK????
        self.prep()

        # Run for 0 fs duration to evaluate system
        my_lmps = self.lammps_system
        my_lmps.run(0)

        # Evaluate potential energy and forces
        pe = my_lmps.eval("pe")
        force_array = []
        for i in range(0, my_lmps.atoms.natoms):
            force_array.append(my_lmps.atoms[i].force)
        return {'potential_energy': pe * u.kcalpermol,
                'forces': force_array * u.kcalpermol / u.angstrom}
    
      
    def prep(self, force=False):
    """
        Drive the construction of the LAMMPS simulation
        This will rebuild this OpenMM simulation if: A) it's not built yet, or B)
        there's a new integrator
    """ 
        
        # If current molecule's velocity is not the same as last recorded velocity, create a new system
        if self._last_velocity == None or self._last_velocity != self.mol.velocities:
            self._create_system()
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

        dataPath = _create_lammps_data(tmpdir, parmedmol)
        print 'Created LAMMPS data'

        pylmp.command("units real")
        pylmp.command("atom_style full")
        pylmp.command("pair_style lj/charmm/coul/long 8.0 10.0 10.0")
        pylmp.command("bond_style harmonic")
        pylmp.command("angle_style harmonic")
        pylmp.command("dihedral_style charmm")
        pylmp.command("improper_style harmonic")
        pylmp.command("kspace_style pppm 0.0001")
        pylmp.command("read_data " + dataPath)

        # securely remove temporary filesystem
        os.remove(dataPath)
        os.umask(saved_umask)
        os.rmdir(tmpdir)
        
        # NOTE: Do we want this?
        pylmp.command("neighbor " + str(self.params.skin) + " bin")
        
        pylmp.command("thermo_style custom step temp pe etotal")
        # NOTE: By default, time step is 10.0
        pylmp.command("thermo 10")
        
        # group hbonds
        hbond_command = _group_hbonds(parmedmol)
        if len(hbond_command) > 0:
            pylmp.command(hbond_command)
            self.group_hbond = True
        else:
            self.group_hbond = False

        # group water
        water_command = _group_water(parmedmol)
        if len(water_command) > 0:
            pylmp.command(water_command)
            self.group_water = True
        else:
            self.group_water = False

        self.lammps_system = pylmp
        
        self._last_velocity = self.mol.velocities.copy() #Keep track of last velocity that was used to create the LAMMPS system 

    
    def _group_hbonds(parmedmol):
    """
        Get indices of non-bond atom types that contain hydrogen bonds
        
        Arguments:
            parmed: parmed struct to iterate through all nonbond types to find hbonds

    """
        hbond_group = ""
        
        if len(parmedmol.LJ_types) <= 0:
            return hbond_group

        for nonbond_name, nonbond_idx in parmedmol.LJ_types.iteritems():
            if(nonbond_name[0] == 'H'):
                hbond_group = hbond_group + " " + str(nonbond_idx)
        
        return "group hbond type" + hbond_group


    def _group_water(parmedmol):
    """
        Get indices of residues that are type 'water'
       
        Arguments:
            parmed: parmed struct to iterate through all residues to find water residues

    """
    
        min_index = sys.maxint;
        max_index = 0;
        for res in parmedmol.residues:
            if res.name == 'WAT':
                if res.idx < min_index:
                    min_index = res.idx
                elif res.idx > max_index:
                    max_index = res.idx
        
        water_group = ""
        if min_index != sys.maxint :
            water_group = "group water molecule <> {0} {1}" .format(min_index, max_index)

        return water_group



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

    
    def _create_lammps_data(tmpdir, parmedmol):
    """
        Create a data. file that contains all necessary information regarding the molecule
        in order to create the appropirate LAMMPS system for it
        
        Arguments:
            tmpdir: temporary directory to store the data file
            parmed: parmed struct to iterate through all nonbond types to find hbonds

    """    
        predictable_filename = 'data.lammps_mol'
        path = os.path.join(tmpdir, predictable_filename)

        with open(path, "w") as lammps_data:   

            # print description
            lammps_data.write("LAMMPS Description\r\n\r\n")
            lammps_data.write("{0} atoms\r\n" .format(len(parmedmol.atoms)))
            lammps_data.write("{0} bonds\r\n" .format(len(parmedmol.bonds)))
            lammps_data.write("{0} angles\r\n" .format(len(parmedmol.angles)))
            lammps_data.write("{0} dihedrals\r\n" .format(len(parmedmol.dihedrals)))
            lammps_data.write("{0} impropers\r\n" .format(len(parmedmol.impropers)))
            
            lammps_data.write("\r\n")
            
            # Get non-bond atom types
            max_idx = 0
            for nonbond_idx in parmedmol.LJ_types.itervalues():
                if(nonbond_idx > max_idx):
                    max_idx = nonbond_idx

            lammps_data.write("{0} atom types\r\n" .format(max_idx))
            lammps_data.write("{0} bond types\r\n" .format(len(parmedmol.bond_types)))
            lammps_data.write("{0} angle types\r\n" .format(len(parmedmol.angle_types)))
            lammps_data.write("{0} dihedral types\r\n" .format(len(parmedmol.dihedral_types)))
            lammps_data.write("{0} improper types\r\n" .format(len(parmedmol.improper_types)))
            
            lammps_data.write("\r\n")

            # calculate lo and hi coordinates for Box size
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

            lammps_data.write("{0} {1} xlo xhi\r\n" .format(xlo, xhi))
            lammps_data.write("{0} {1} ylo yhi\r\n" .format(ylo, yhi))
            lammps_data.write("{0} {1} zlo zhi\r\n" .format(zlo, zhi))
            
            lammps_data.write("\r\n")


            # print masses -- Number of atom types
            lammps_data.write("Masses\r\n\r\n")
            for i in range(1, max_idx+1):
                for atom in parmedmol.atoms:
                    if atom.nb_idx == i:
                        lammps_data.write("{0} {1}\r\n" .format(i, atom.mass))
                        break
            
            lammps_data.write("\r\n")
            

            # Pair Coeffs -- Number of atom types
            lammps_data.write("Pair Coeffs\r\n\r\n")
            for i in range(1, max_idx+1):
                for atom in parmedmol.atoms:
                    if atom.nb_idx == i:
                        lammps_data.write("{0} {1} {2} {3} {4}\r\n" .format(i, atom.epsilon, atom.sigma, 
                                                                            atom.epsilon_14, atom.sigma_14))
                        break
            
            lammps_data.write("\r\n")
            

            # Bond Coeffs - Harmonic
            if parmedmol.bond_types:
                lammps_data.write("Bond Coeffs\r\n\r\n")
                for bt in parmedmol.bond_types:
                    lammps_data.write("{0} {1} {2}\r\n" .format(bt.idx+1, bt.k, bt.req))
                
                lammps_data.write("\r\n")


            # Angle Coeffs - Harmonic
            if parmedmol.angle_types:
                lammps_data.write("Angle Coeffs\r\n\r\n")
                for at in parmedmol.angle_types:
                    lammps_data.write("{0} {1} {2}\r\n" .format(at.idx+1, at.k, at.theteq))
               
                lammps_data.write("\r\n")


            # Dihedral Coeffs - Charmm
            if parmedmol.dihedral_types:
                lammps_data.write("Dihedral Coeffs\r\n\r\n")
                for dt in parmedmol.dihedral_types:
                    lammps_data.write("{0} {1} {2} {3} {4}\r\n" .format(dt.idx+1, dt.phi_k, dt.per, 
                        (180 if dt.phase >= 180 else 0), 
                        (1.0 if dt.scnb == 1.0 else 0.0)))
                
                lammps_data.write("\r\n")


            # Improper Coeffs - Harmonic
            if parmedmol.improper_types:
                lammps_data.write("Improper Coeffs\r\n\r\n")
                for it in parmedmol.improper_types:
                    print
                    lammps_data.write("{0} {1} {2}\r\n" .format(it.idx+1, it.psi_k, it.psi_eq))
                
                lammps_data.write("\r\n")


            # print atoms
            lammps_data.write("Atoms\r\n\r\n")
            for atom in parmedmol.atoms:
                lammps_data.write("{0} {1} {2} {3} " .format(atom.idx+1, atom.residue.idx, atom.nb_idx, atom.charge))
                if not self.mol.positions.any():        
                    lammps_data.write("{0} {1} {2} " .format(0.00, 0.00, 0.00))        
                else:
                    pos = self.mol.positions[atom.idx]
                    lammps_data.write("{0} {1} {2} " .format(pos[0].value_in(u.angstrom), pos[1].value_in(u.angstrom), 
                        pos[2].value_in(u.angstrom)))
                lammps_data.write("{0} {1} {2}\r\n" .format(0, 0, 0)) # nx, ny, nz not defined    

            lammps_data.write("\r\n")


            # print velocities
            lammps_data.write("Velocities\r\n\r\n")
            for atom in self.mol.atoms:
                xvel = atom.vx.value_in(u.angstrom/u.fs)
                yvel = atom.vy.value_in(u.angstrom/u.fs)
                zvel = atom.vz.value_in(u.angstrom/u.fs)
                lammps_data.write("{0} {1} {2} {3}\r\n" .format(atom.index+1, xvel, yvel, zvel))
                
            lammps_data.write("\r\n")


            # print bonds
            lammps_data.write("Bonds\r\n\r\n")
            for idx, bond in enumerate(parmedmol.bonds):
                lammps_data.write("{0} {1} {2} {3}\r\n" .format(idx+1, bond.type.idx+1, 
                  bond.atom1.idx+1, bond.atom2.idx+1))
                
            lammps_data.write("\r\n")


            # print angles
            lammps_data.write("Angles\r\n\r\n")
            for idx, angle in enumerate(parmedmol.angles):
                lammps_data.write("{0} {1} {2} {3} {4}\r\n" .format(idx+1, angle.type.idx+1,
                  angle.atom1.idx+1, 
                  angle.atom2.idx+1, angle.atom3.idx+1))
                
            lammps_data.write("\r\n")


            # print dihedrals
            if parmedmol.dihedrals:
                lammps_data.write("Dihedrals\r\n\r\n")
                for idx, di in enumerate(parmedmol.dihedrals):
                    lammps_data.write("{0} {1} {2} {3} {4} {5}\r\n" .format(idx+1, di.type.idx+1, 
                      di.atom1.idx+1, di.atom2.idx+1, 
                      di.atom3.idx+1, di.atom4.idx+1))
                    
                lammps_data.write("\r\n")


            # print impropers
            if parmedmol.impropers:
                lammps_data.write("Impropers\r\n\r\n")
                for im in parmedmol.impropers:
                    lammps_data.write("{0} {1} {2} {3} {4} {5}\r\n" .format(im.idx+1, im.type.idx+1, 
                      im.atom1.idx+1, im.atom2.idx+1, 
                      im.atom3.idx+1, im.atom4.idx+1))
                
                lammps_data.write("\r\n")

        
        return path


