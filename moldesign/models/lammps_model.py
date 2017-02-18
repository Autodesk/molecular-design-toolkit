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
    """Creates an OpenMM "context" to drive energy calculations.
    Note that, while a dummy integrator is assigned, a different context will
    be created for any MD calculations.

    :ivar sim: openmm simulation object
    :type sim: simtk.openmm.app.Simulation
    """
    # NEWFEATURE: need to set/get platform (and properties, e.g. number of threads)
    DEFAULT_PROPERTIES = ['potential_energy', 'forces']

    def __init__(self, **kwargs):
        super(LAMMPSPotential, self).__init__(**kwargs)
        self._prepped = False # is the model prepped?
        
        # Other classes have access to these attributes
        self.lammps_system = None # Basic LAMMPS system for this molecule
        self.hbond_group = None # Hydrogen mass 
        self.unit_system = None
        self.export_to_web = None
        
    # calculate potential energy and force
    def calculate(self, requests):
        
        # Recreate system if molecule changed
        results = self.evaluate_molecule()
        return results

    def evaluate_molecule(self):
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

    def prep(self):
        """
        Drive the construction of the LAMMPS simulation
        This will rebuild this OpenMM simulation if: A) it's not built yet, or B)
        there's a new integrator
        """ 

        if not self._prepped or self.lammps_system is None:
            self.lammps_system = self._create_system()
            
            self.unit_system = u.UnitSystem(length=u.angstrom, force=u.kcalpermol/u.fs, energy=u.kcalpermol, time=u.fs,
                mass=u.gpermol, charge=u.electron_charge)

        self._prepped = True

    #################################################
    # "Private" methods for managing LAMMPS are below
    
    def _create_system(self):

        """
        Create a LAMMPS system. Use MDT molecule object to construct a LAMMPS system

        """

        # Ensure force fields are assigned 
        if 'amber_params' not in self.mol.ff:
            raise NotImplementedError('Assign forcefield parameters to the system')

        # Store the PRMTOP file if ff are assigned and use it to create lammps molecule data

        # initialize lammps object
        lmp = lammps()
        pylmp = PyLammps(ptr=lmp)

        # create temporary file system
        tmpdir = tempfile.mkdtemp()
        saved_umask = os.umask(0077)
        fromatted_data = self.format_lammps_data()

        # Create temperary data file
        predictable_filename = 'data.lammps_mol'
        data_path = os.path.join(tmpdir, predictable_filename)

        with open(data_path, "w") as lammps_data:
            lammps_data.write(fromatted_data)

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

        pylmp.command("neighbor    2.0 bin")
        # pylmp.command("neigh_modify    delay 5")

        # securely remove temporary filesystem
        os.remove(data_path)
        os.umask(saved_umask)
        os.rmdir(tmpdir)

        return pylmp

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
    
    def format_lammps_data(self):
        """
        Parse molecule and force fields info to generate a formatted data

        """    
        # Create parmed object
        self.mol.ff.amber_params.prmtop.put('/tmp/tmp.prmtop')
        parmedmol = med.load_file('/tmp/tmp.prmtop')

         # Get non-bond atom types
        parmedmol.fill_LJ()
        max_idx = 0
        for nonbond_idx in parmedmol.LJ_types.itervalues():
            if nonbond_idx > max_idx:
                max_idx = nonbond_idx

        datalines =  ""
        atom_types = ""

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

        datalines += "{0} {1} xlo xhi\r\n".format(-200, 200)
        datalines += "{0} {1} ylo yhi\r\n".format(-200, 200)
        datalines += "{0} {1} zlo zhi\r\n".format(-200, 200)

        datalines += "\r\n"

        # Masses - atom types
        datalines += "Masses\r\n\r\n"
        for i in range(1, max_idx+1):
            for atom in parmedmol.atoms:
                if atom.atomic_number == 1 and self.hbond_group is None:
                    self.hbond_group = atom.mass
                if atom.nb_idx == i:
                    atom_types += "{0} ".format(self.mol.atoms[atom.idx].symbol)
                    datalines += "{0} {1}\r\n".format(i, atom.mass)
                    break

        atom_types += "\r\n"
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
        self.export_to_web = atom_types + datalines
        # print self.export_to_web
        return datalines

@exports
class LAMMPSInteractive(EnergyModelBase):
    NAME_AFFECTED_ATOMS = "affected_atom"
    NAME_NEARBY_ATOMS = "nearby_atoms"
    NUM_AFFECTED_RADIUS = 5

    PARAMETERS = [
        mdt.parameters.Parameter(
            'affected_atom', 'Atom affected by user\'s interaction', type=mdt.molecules.Atom, default=None)]

    def __init__(self, **kwargs):
        super(LAMMPSInteractive, self).__init__(**kwargs)
        self._prepped = False   # is the model prepped?
        self._potential_model = None
        self._calculated_properties = None

        self._affected_atom = None  # atom that the user is directly interacting with
        self._force_vector = None   # force that acts on the molecule
        self._initial_position = None   # initial position of atom that the user is interacting with
        self._saved_positions = []  # Saved positions of the affected atom and its nearby atoms
        self._desired_position = None   # desired position of the affected atom
        self._stretch_factor = None     # stretch factor of the affected atom's neighbor atoms due to indirect force

        self.lammps_system = None   # LAMMPS object for this molecule
        self.hbond_group = None
        self.unit_system = None

    def calculate(self, requests):
        '''
            After assigning this energy model to the molecule, call this function as soon as possible, so
            the molecule is ready for the interactive simulation
        '''
        
        # Recreate system if molecule changed
        if self._calculated_properties is None:
            self.prep()
            self._calculated_properties = self._potential_model.evaluate_molecule()

        return self._calculated_properties

    def render_molecule(self, traj=None):
        '''
            For interactive purposes, reset the atom of interest to desired position based on force applied by user
            Call this function after each round of simulation
        '''
        if self._desired_position is None:
            raise ValueError("Can't call reset position before running any interactive simulation!")

        # Save positions after simulation
        for idx in range(0, self.NUM_AFFECTED_RADIUS*2+1):
            nearby_atom_pos = self.mol.positions[self._affected_atom.index - self.NUM_AFFECTED_RADIUS + idx].copy()
            if traj is not None:
                nearby_atom_pos = traj.positions[self._affected_atom.index - self.NUM_AFFECTED_RADIUS + idx].copy()

            self._saved_positions.insert(idx, nearby_atom_pos)

        # render appropriate positions to affected atom + nearby atoms
        self.mol.positions[self._affected_atom.index] = self._desired_position
        if traj is not None:
            traj.positions[self._affected_atom.index] = self._desired_position

        for i in xrange(self.NUM_AFFECTED_RADIUS):
            self.mol.positions[self._affected_atom.index - (i + 1)] += self._stretch_factor.copy()
            self.mol.positions[self._affected_atom.index + (i + 1)] += self._stretch_factor.copy()

            if traj is not None:
                traj.positions[self._affected_atom.index - (i + 1)] += self._stretch_factor.copy()
                traj.positions[self._affected_atom.index + (i + 1)] += self._stretch_factor.copy()

    def finish_interaction(self):
        '''
            When the user is done interacting with the molecule (i.e. lets go of the mouse), 
            reset the atom's position to chemicaly-correct position

            Prior to next interaction, reset the model
        '''
        if self._saved_positions is None:
            raise ValueError("Can't finish interaction if haven't run any!")

        # Save positions after simulation
        for idx in range(0, self.NUM_AFFECTED_RADIUS*2+1):
            saved_pos = self._saved_positions[idx];
            self.mol.positions[self._affected_atom.index - self.NUM_AFFECTED_RADIUS + idx] = saved_pos

    def apply_force(self, new_pos):
        '''
            :param new_pos: amount by which we want to change the selected atom's position
            Call this function each time user applies a different force (i.e. move the mouse cursor somewhere else)
        '''
        pos_vector = new_pos.value_in(u.angstrom)
        pointer_vector = [10 * p / abs(p) if p != 0.0 else p for p in pos_vector]
        vel_vector = [p / 65 * 1.0 for p in pos_vector]
        mass = self._affected_atom.mass.value_in(u.amu)
        self._force_vector = [v * abs(v) * mass / 2 for v in vel_vector]
        self._desired_position = self._initial_position + (pointer_vector * u.angstrom)
        self._stretch_factor = [p / 2 for p in pointer_vector] * u.angstrom

    def prep(self):
        """
        Drive the construction of the LAMMPS simulation
        This will rebuild this OpenMM simulation if: A) it's not built yet, or B)
        there's a new integrator
        """

        if self.params.affected_atom is None:
            raise ValueError("Must initialize energy model with an atom for interaction")

        if not self._prepped or self._potential_model is None or self.lammps_system is None:
            self._potential_model = mdt.models.LAMMPSPotential()
            self._potential_model.mol = self.mol
            self._potential_model.prep()
            self._calculated_properties = None
            self.lammps_system = self._potential_model.lammps_system

            # Group affected atom
            lmps = self.lammps_system
            self._affected_atom = self.params.affected_atom
            self._initial_position = self.mol.positions[self._affected_atom.index].copy()
            group_atom = "group {0} id {1}:{2}".format(self.NAME_AFFECTED_ATOMS,
                                                       self._affected_atom.index - self.NUM_AFFECTED_RADIUS + 1,
                                                       self._affected_atom.index + self.NUM_AFFECTED_RADIUS + 1)
            lmps.command(group_atom.rstrip())

            # Get settings from lammps model
            self.unit_system = self._potential_model.unit_system
            self.hbond_group = self._potential_model.hbond_group

        # Apply user force fix
        if self._force_vector is not None:
            lmps = self.lammps_system
            force_vector = self._force_vector
            print force_vector
            # lmps.command("fix apply_user_force {0} addforce {1} {2} {3}"
            #              .format(self.NAME_AFFECTED_ATOMS, force_vector[0], force_vector[1], force_vector[2]))

            lmps.command("fix apply_user_force all addforce {0} {1} {2}".format(force_vector[0],
                                                                                force_vector[1],
                                                                                force_vector[2]))

        # Dump positions only once throughout the simulation

        # Reset to chemically stable position before resuming interactive simulation
        if self._saved_positions:
            for idx in range(0, self.NUM_AFFECTED_RADIUS*2+1):
                saved_pos = self._saved_positions[idx]
                self.mol.positions[self._affected_atom.index - self.NUM_AFFECTED_RADIUS + idx] = saved_pos

        self._prepped = True

    def get_data_for_web(self):
        self.prep()
        if self._potential_model.export_to_web is not None:
            return self._potential_model.export_to_web





