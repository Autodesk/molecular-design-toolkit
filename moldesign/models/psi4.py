from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library

standard_library.install_aliases()

# Copyright 2017 Autodesk Inc.
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

from .. import units as u
from .. import orbitals
from .. import molecules
from ..compute import runsremotely
from ..interfaces import psi4_interface
from .base import QMBase

# On June 15, JSOB added " import numpy as np
import numpy as np

# On June 13, JSOB added " ,'mp2':'mp2' "

THEORIES = {'rhf':'hf', 'mp2':'mp2', 'uhf':'hf', 'dft':'b3lyp', 'casscf':'casscf'}

BASES = {}

# Names for one-electron properties as mdt_name : psi4_name

ONE_ELECTRON_PROPERTIES = { 'dipole':'DIPOLE',
                            'quadrupole':'QUADRUPOLE',
                            'electrostatic_potential_at_nuclei' : 'ESP_AT_NUCLEI' ,
                            'molecular_orbital_extents' : 'MO_EXTENTS' ,
                            'mulliken_charges': 'MULLIKEN_CHARGES',
                            'lowdin_charges' : 'LOWDIN_CHARGES',
                            'wiberg_lowdin_indices': 'WIBERG_LOWDIN_INDICES',
                            'mayer_indices':'MAYER_INDICES' ,
                            'natural_orbital_occupations':'NO_OCCUPATIONS' }

# NATIVE_OPTIONS added by JSOB
NATIVE_OPTIONS = {'mdt_frozen_core':'freeze_core'} # To be edited as MDT develops


# returns mdt trajectory
def capture_history(name_string, mdt_molecule, psi4_history, model):
    import psi4
    my_traj=molecules.Trajectory(mol=mdt_molecule)
    my_traj.__setattr__('molecules' , [])
        
    for step in range(len(psi4_history['energy'])):
        coords = np.asarray(psi4_history['coordinates'][step])
        atoms = []
        for atom_count in range(len(mdt_molecule.atoms)):
            this_atom = molecules.Atom(name    = mdt_molecule.atoms[atom_count].name,
                             element = mdt_molecule.atoms[atom_count].element
                            )
            this_atom.x = coords[atom_count][0] * u.angstrom
            this_atom.y = coords[atom_count][1] * u.angstrom
            this_atom.z = coords[atom_count][2] * u.angstrom
            atoms.append(this_atom)
        
        step_molecule=molecules.Molecule( atoms, name = name_string + str(step))
        step_molecule.set_energy_model(model)
        step_molecule.properties['potential_energy'] = psi4_history['energy'][step]
        my_traj.frames.append(step_molecule.properties)
        my_traj.molecules.append(step_molecule)
    my_traj.__setattr__('potential_energy', psi4_history['energy'])

    print(my_traj)
    
    return my_traj


class Psi4Potential(QMBase):
    def prep(self):
        return True

    @runsremotely(enable=psi4_interface.force_remote)
    def calculate(self, requests):
        import psi4
        
        print(requests)
        
        request_options = {'vibrational_frequencies' : psi4.freq, 'freq' : psi4.freq, 'frequency' : psi4.freq,'frequencies' : psi4.freq,'optimized_geometry' : psi4.opt,'opt' : psi4.opt,'optimize' : psi4.opt,'potential_energy' : psi4.energy,'energy' : psi4.energy,'electronic_gradient' : psi4.gradient,'hessian' : psi4.hessian }

        try:
            assert self.params.dimer_partner
        except AttributeError:
            psi4_molecule = psi4_interface.mdt_to_psi4(self.mol)
        else:
            psi4_molecule = psi4_interface.mdt_to_psi4_dimer(self.mol, self.params.dimer_partner)


        psi4_molecule = psi4_interface.mdt_to_psi4(self.mol)
        
        self.name = psi4_interface.mol_name(self.mol)
        
        self._psi4_set_up(requests)
        
        if 'memory' in self.params:
            self._set_memory()

        his=False
        psi4_results =self._run_psi4_calc(requests, psi4_molecule, his)
        
        props=self._get_properties(requests, psi4_results)

        psi4.core.clean()
        psi4.core.clean_options()
        return props

    @runsremotely(enable=psi4_interface.force_remote)
    def minimize(self, requests):
        import psi4
        self._psi4_set_up(requests)
        
        
        if 'memory' in self.params:
            self._set_memory()
        
        try:
            assert self.params.dimer_partner
        except AttributeError:
            psi4_molecule = psi4_interface.mdt_to_psi4(self.mol)
        else:
            psi4_molecule = psi4_interface.mdt_to_psi4_dimer(self.mol, self.params.dimer_partner)

        self.name = psi4_interface.mol_name(self.mol)
        
        self._psi4_set_up(requests)
        
        his=True
        psi4_results = self._run_psi4_calc(requests, psi4_molecule, his)
        
        prop_output=[psi4_results[0],psi4_results[1]]
        self.mol.props = self._get_properties(requests, prop_output )
        
        print(requests)
        
        trajectory=capture_history(self.name, self.mol, psi4_results[2], self)
        
        psi4.core.clean()
        psi4.core.clean_options()
        return trajectory

    # JSOB added June 19
    def _set_memory(self):
        import psi4
        psi4.set_memory(self.params.memory)
    
    def _psi4_set_up(self, requests):
        import psi4

        request_options = {'vibrational_frequencies' : psi4.freq, 'freq' : psi4.freq, 'frequency' : psi4.freq,'frequencies' : psi4.freq,'optimized_geometry' : psi4.opt,'opt' : psi4.opt,'optimize' : psi4.opt,'potential_energy' : psi4.energy,'energy' : psi4.energy,'electronic_gradient' : psi4.gradient,'hessian' : psi4.hessian }
        
        try:
            assert requests
        except AssertionError:
            pass
        else:
            for req in requests:
                if req in request_options.keys():
                    if request_options[req] == psi4.freq:
                        self.params.options.update({'':''})

        self.params.theory=self._gettheorystring()
        try:
            assert self.params.options
        except AttributeError:
            pass
        else:
            psi4.set_options(self._getoptionsstring())
        
        try:
            assert self.params.output
        except AttributeError:
            pass
        else:
            psi4.core.set_output_file(self.params.output, False)
        
        try:
            assert self.params.basis
        except AssertionError:
            pass
        else:
            psi4.set_options({'basis': BASES.get(self.params.basis, self.params.basis)})

        try:
            assert self.params.dertype
        except AttributeError:
            self.params.dertype = 'gradient'
        else:
            pass
            
    
    @staticmethod
    def _get_runmethod(requests):
        import psi4
        request_options = {'vibrational_frequencies' : psi4.freq, 'freq' : psi4.freq, 'frequency' : psi4.freq,'frequencies' : psi4.freq ,'optimized_geometry' : psi4.opt,'opt' : psi4.opt,'optimize' : psi4.opt,'potential_energy' : psi4.energy,'energy' : psi4.energy,'electronic_gradient' : psi4.gradient,'hessian' : psi4.hessian, 'gradient':psi4.gradient }
        
        for quantity in requests:
            if quantity == 'potential_energy':
                continue
            if quantity in request_options.keys():
                return request_options[quantity]
            else:
                pass


        print('Returning potential_energy from psi4.energy by default')
        return psi4.energy

    def _gettheorystring(self):
        import psi4
        
        self.params.theory=self.params.theory.lower()
        
        try:
            assert self.params.options
        except AttributeError:
            if self.params.theory in ('uhf', 'UHF'):
                self.params.options={'reference': 'uhf'}
            else:
                pass
        else:
            if self.params.theory in ('uhf', 'UHF'):
                self.params.options.append({'reference':'uhf'})
        
        return THEORIES.get(self.params.theory, self.params.theory)

    def _getoptionsstring(self, requests):
        for mdt_opt in self.params.options:
            if mdt_opt in NATIVE_OPTIONS:
                self.params.options[NATIVE_OPTIONS[mdt_opt]]=self.params.options[mdt_opt]
                del self.params.options[mdt_opt]
        return self.params.options
    
    def _run_psi4_calc(self, requests, psi4_molecule, his):
        import psi4
        print(requests)
        psi4_output = self._get_runmethod(requests)(self.params.theory, molecule=psi4_molecule, return_wfn=True, return_history=his )
        psi4_molecule.update_geometry()
        self.mol.charge=psi4_molecule.molecular_charge()
        self.mol.multiplicity=psi4_molecule.multiplicity()
        return psi4_output
    
    #psi4_output should be a two or three membered list: [psi4_energy, psi4_wavefunction] or [ psi4_energy, psi4_wavefunction, psi4_history ]
    def _get_properties(self,requests,psi4_output):
        import psi4
        props={}
        
        props['potential_energy']=psi4_output[0]
    
    #DO NOT DELETE LINE 206 WE NEED IT FOR TESTING
        props['psi4_output'] = psi4_output
        psi4_vars_title=self.name
        
        psi4.oeprop(psi4_output[1], 'DIPOLE',title=psi4_vars_title)
        psi4_variables=psi4.core.get_variables()
        
        #DEFAULT_PROPERTIES evaluation
        props['nuclear_repulsion'] = psi4_variables['NUCLEAR REPULSION ENERGY']
        props['dipole_moment']     = np.linalg.norm([ psi4_variables[psi4_vars_title +' DIPOLE '+x]  for x in ('X', 'Y', 'X') ]) * u.debye
        
        props.update(self._get_one_e_props(requests, psi4_output[1],psi4_vars_title ))
        
        
        
        #props.update(self._get_ao_basis_functions(self, requests, psi4_output[1], psi4_vars_title))
        #props['wfn'] = self._build_mdt_wavefunction(psi4_output[1])
        
        #Capture History
        return props
########################################################
# DO NOT END _get_properties ABOVE OR BELOW THIS LINE
########################################################

    def _get_one_e_props(self, requests, psi4_wavefunction, psi4_vars_title):
        oeprops = {}
        oeprop_funcs = {'dipole':self.get_dipole,
'quadrupole':self.get_quadrupole,
'electrostatic_potential_at_nuclei':self.get_esp_at_nuclei,
'molecular_orbital_extents':self.get_mo_extents,
'mulliken_charges':self.get_mulliken_charges,
'lowdin_charges':self.get_lowdin_charges,
'wiberg_lowdin_indices':self.get_wib_indices,
'mayer_indices':self.get_mbi_indices,
'natural_orbital_occupations':self.get_no_occupations}
    
        for req in requests:
            if req in oeprop_funcs.keys():
                oeprops.update({req : oeprop_funcs[req]( psi4_wavefunction, psi4_vars_title) })
        
        return oeprops

    def get_dipole(self, psi4_wavefunction,psi4_vars_title):
        import psi4
        psi4.oeprop(psi4_wavefunction, 'DIPOLE', title =psi4_vars_title)
        psi4_variables = psi4.core.get_variables()
        dip = [ psi4_variables[psi4_vars_title +' DIPOLE '+x]  for x in ('X', 'Y', 'X') ]
        return dip
    
    def get_quadrupole(self, psi4_wavefunction,psi4_vars_title):
        import psi4
        psi4.oeprop(psi4_wavefunction, 'QUADRUPOLE', title = psi4_vars_title)
        psi4_variables = psi4.core.get_variables()
        qup = [ psi4_variables[psi4_vars_title+' QUADRUPOLE ' + x] for x in ('XX','XY','XZ','YY','YY','ZZ') ]
        return qup
    
    def get_esp_at_nuclei(self, psi4_wavefunction,psi4_vars_title):
        import psi4
        psi4.oeprop(psi4_wavefunction,'ESP_AT_NUCLEI', title = psi4_vars_title )
        psi4_variables = psi4.core.get_variables()
        esp = []
        
        for i in range(1,(self.mol.num_atoms+1)):
            esp.append(psi4_variables['ESP AT CENTER '+str(i)])
        return esp
    
    def get_mo_extents(self, psi4_wavefunction,psi4_vars_title):
        import psi4
        
        return moe
    
    def get_mulliken_charges(self, psi4_wavefunction,psi4_vars_title):
        import psi4
        import numpy as np
        psi4.oeprop(psi4_wavefunction, 'MULLIKEN_CHARGES')
        nuq = np.asarray(psi4_wavefunction.atomic_point_charges())
        return nuq
    
    def get_lowdin_charges(self, psi4_wavefunction,psi4_vars_title):
        import psi4
        import numpy as np
        psi4.oeprop(psi4_wavefunction, 'LOWDIN_CHARGES')
        loq = np.asarray(psi4_wavefunction.atomic_point_charges())
        return loq
    
    def get_wib_indices(self,psi4_wavefunction,psi4_vars_title):
        import psi4
        psi4.oeprop(psi4_wavefunction,'WIBERG_LOWDIN_INDICES')
        wbi = np.asarray(psi4_wavefunction.get_array('WIBERG_LOWDIN_INDICES'))
        return wbi
    
    def get_mbi_indices(self,psi4_wavefunction,psi4_vars_title):
        import psi4
        psi4.oeprop(psi4_wavefunction, 'MAYER_INDICES')
        mbi = np.asarray(psi4_wavefunction.get_array('MAYER_INDICES'))
        return mbi
    
    def get_no_occupations(self,psi4_wavefunction,psi4_vars_title):
        import psi4
        psi4.oeprop(psi4_wavefunction, 'NO_OCCUPATIONS' )
        noo = psi4_wavefunction.no_occupations()
        print(noo)
        return noo

    def _build_mdt_wavefunction(self, psi4_wavefunction):
        import psi4
        
        bfs = []
                              
        specd_basis = orbitals.basis.BasisSet(mol=self.mol,
                                              orbitals = self._get_ao_basis_functions(psi4_wavefunction),
                                              name = psi4_wavefunction.basisset().name(),
                                              h1e = np.asarray(psi4_wavefunction.H()),
                                              overlaps = np.asarray(psi4_wavefunction.S)
                                              )
        
        returned_wfn = orbitals.wfn.ElectronicWfn(mol = self.mol)
    
        return returned_wfn

    def _get_ao_basis_functions(self, psi4_wavefunction):
        import psi4
        from psi4 import qcdb

        psi4_molecule = psi4_interface.mdt_to_psi4(self.mol)

        mymol = qcdb.Molecule(psi4_molecule.create_psi4_string_from_molecule())
        
        print('[1]    <<<  uniform cc-pVDZ  >>>')
        wert = qcdb.BasisSet.pyconstruct(mymol, 'BASIS', 'cc-pvdz')[0]
        print(wert.print_detail())
        psi4.compare_integers(38, wert.nbf(), 'nbf()')
        psi4.compare_integers(40, wert.nao(), 'nao()')
        psi4.compare_strings('c2v', mymol.schoenflies_symbol(), 'symm')
        
        bfs = []
        #      pmol = self.pyscfmol
        
        orblabels = iter(pmol.spheric_labels())
        
        # Rewrite to loop over qcdb shells making sure to
        for ishell in range(wert.nbf()):  # loop over shells (n,l)
            atom = my_mol.atoms[pmol.bas_atom(ishell)]
            angular = pmol.bas_angular(ishell)
            num_momentum_states = angular*2 + 1
            exps = pmol.bas_exp(ishell)
            num_contractions = pmol.bas_nctr(ishell)
            coeffs = pmol.bas_ctr_coeff(ishell)
            
            for ictr in range(num_contractions):  # loop over contractions in shell
                for ibas in range(num_momentum_states):  # loop over angular states in shell
                    label = next(orblabels) # 2 parts orientation & quantum number
                    sphere_label = label[3]
                    l, m = SPHERICAL_NAMES[sphere_label]
                    assert l == angular
                    # TODO: This is not really the principal quantum number
                    n = int(''.join(x for x in label[2] if x.isdigit()))
                    
                    primitives = [orbitals.SphericalGaussian(atom.position.copy(),
                                                             exp, n, l, m,
                                                             coeff=coeff[ictr])
                                  for exp, coeff in zip(exps, coeffs)]
                    bfs.append(orbitals.AtomicBasisFunction(atom, n=n, l=angular, m=m,
                                                                              primitives=primitives))
        
        return bfs


#        for atoms in self.mol.atoms:
#    bfs.append(orbitals.AtomicBasisFunction())














