import psi4
import moldesign as mdt
#This module tests the consistency between calculations submitted by the
#molecular design toolkit and those sumbitted "directly" to psi4. It only
#compares psi4 values and does not compare MDT values to Psi4 values.

class Psi4_to_Psi4_test(object):
    
    def __init__(self, mem = '500 MB', **kwargs):
        self.mem = mem
        self.specs = {}
        self.specs.update(kwargs)
    

    def energy_to_energy(self):
        psi4.set_memory(self.mem)
        
        try:
            assert self.specs['psi4_options']
        except KeyError:
            pass
        else:
            psi4.set_options(self.specs['psi4_options'])
        
        print("Running direct psi4 calculation")
        energy_from_psi4 , wfn_from_psi4 = psi4.energy(self.specs['method_name']+'/'+self.specs['basis_name'], molecule=self.specs['psi4_mol'], return_wfn=True)
        
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()
        
        my_mdt_mol = mdt.interfaces.psi4_interface.psi4_to_mdt(self.specs['psi4_mol'])

        try:
            assert self.specs['psi4_options']
        except KeyError:
            my_mdt_mol.set_energy_model(mdt.models.psi4.Psi4Potential,
                                        basis=self.specs['basis_name'],
                                        memory=self.mem,
                                        theory=self.specs['method_name'])
        else:
            my_mdt_mol.set_energy_model(mdt.models.psi4.Psi4Potential,
                                        basis=self.specs['basis_name'],
                                        memory=self.mem,
                                        options=self.specs['psi4_options'],
                                        theory=self.specs['method_name'])

        print("Running psi4 calculation using MDT")
        energy_from_mdt , wfn_from_mdt = my_mdt_mol.calculate()['psi4_output']
    
        assert psi4.compare_values(energy_from_psi4 , energy_from_mdt, 6 , "energy_comparison")

        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()
    
    def freq_to_freq(self):
        
        psi4.set_memory(self.mem)
        self.specs['psi4_mol'].update_geometry()
        try:
            assert self.specs['psi4_options']
        except KeyError:
            pass
        else:
            psi4.set_options(self.specs['psi4_options'])
        
        print("Running direct psi4 calculation")
        energy_from_psi4 , wfn_from_psi4 = psi4.freq(self.specs['method_name']+'/'+self.specs['basis_name'], molecule=self.specs['psi4_mol'], return_wfn=True)
        freqs_from_psi4 = [ wfn_from_psi4.frequencies().get(i) for i in range(wfn_from_psi4.nirrep())]

        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()
        
        my_mdt_mol = mdt.interfaces.psi4_interface.psi4_to_mdt(self.specs['psi4_mol'])

        try:
            assert self.specs['psi4_options']
        except KeyError:
            my_mdt_mol.set_energy_model(mdt.models.psi4.Psi4Potential,
                                        basis=self.specs['basis_name'],
                                        memory=self.mem,
                                        theory=self.specs['method_name'])
        else:
            my_mdt_mol.set_energy_model(mdt.models.psi4.Psi4Potential,
                                        basis=self.specs['basis_name'],
                                        memory=self.mem,
                                        options=self.specs['psi4_options'],
                                        theory=self.specs['method_name'])
        
        print("Running psi4 calculation using MDT")
        energy_from_mdt, wfn_from_mdt = my_mdt_mol.calculate(requests = ['freq'])['psi4_output']
        freqs_from_mdt = [wfn_from_mdt.frequencies().get(i) for i in range(wfn_from_mdt.nirrep())]

        assert wfn_from_mdt.nirrep() == wfn_from_psi4.nirrep()
        
        for i in range(wfn_from_psi4.nirrep()):
            assert psi4.compare_values(freqs_from_psi4[i] , freqs_from_mdt[i], 1, "Frequency comparison")
        
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()

    def hessian_to_hessian(self):
        psi4.set_memory(self.mem)
        self.specs['psi4_mol'].update_geometry()
        try:
            assert self.specs['psi4_options']
        except KeyError:
            pass
        else:
            psi4.set_options(self.specs['psi4_options'])
        
        print("Running direct psi4 calculation")
        hessian_from_psi4 , wfn_from_psi4 = psi4.hessian(self.specs['method_name']+'/'+self.specs['basis_name'], molecule=self.specs['psi4_mol'], return_wfn=True)
        
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()

        my_mdt_mol = mdt.interfaces.psi4_interface.psi4_to_mdt(self.specs['psi4_mol'])
        
        try:
            assert self.specs['psi4_options']
        except KeyError:
            my_mdt_mol.set_energy_model(mdt.models.psi4.Psi4Potential,
                                        basis=self.specs['basis_name'],
                                        memory=self.mem,
                                        theory=self.specs['method_name'])
        else:
            my_mdt_mol.set_energy_model(mdt.models.psi4.Psi4Potential,
                                        basis=self.specs['basis_name'],
                                        memory=self.mem,
                                        options=self.specs['psi4_options'],
                                        theory=self.specs['method_name'])
                                    
        print("Running psi4 calculation using MDT")
        hessian_from_mdt, wfn_from_mdt = my_mdt_mol.calculate(requests = ['hessian'])['psi4_output']
        psi4.compare_matrices( hessian_from_psi4, hessian_from_mdt , 6, self.specs['method_name']+"/"+self.specs['basis_name']+" Hessian" )

        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()

    def gradient_to_gradient(self):
        import numpy as np
        psi4.set_memory(self.mem)
        self.specs['psi4_mol'].update_geometry()
        try:
            assert self.specs['psi4_options']
        except KeyError:
            pass
        else:
            psi4.set_options(self.specs['psi4_options'])
        
        print("Running direct psi4 calculation")
        gradient_from_psi4 , wfn_from_psi4 = psi4.gradient(self.specs['method_name']+'/'+self.specs['basis_name'], molecule=self.specs['psi4_mol'], return_wfn=True)

        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()
        
        my_mdt_mol = mdt.interfaces.psi4_interface.psi4_to_mdt(self.specs['psi4_mol'])
        
        try:
            assert self.specs['psi4_options']
        except KeyError:
            my_mdt_mol.set_energy_model(mdt.models.psi4.Psi4Potential,
                                        basis=self.specs['basis_name'],
                                        memory=self.mem,
                                        theory=self.specs['method_name'])
        else:
            my_mdt_mol.set_energy_model(mdt.models.psi4.Psi4Potential,
                                        basis=self.specs['basis_name'],
                                        memory=self.mem,
                                        options=self.specs['psi4_options'],
                                        theory=self.specs['method_name'])
        
        print("Running psi4 calculation using MDT")
        gradient_from_mdt, wfn_from_mdt = my_mdt_mol.calculate(requests = ['gradient'])['psi4_output']
        psi4.compare_matrices( gradient_from_psi4, gradient_from_mdt , 6, self.specs['method_name']+"/"+self.specs['basis_name']+" Hessian" )

        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()

    def optimize_to_optimize(self):
        psi4.set_memory(self.mem)
        self.specs['psi4_mol'].update_geometry()
        my_mdt_mol = mdt.interfaces.psi4_interface.psi4_to_mdt(self.specs['psi4_mol'])
        try:
            assert self.specs['psi4_options']
        except KeyError:
            pass
        else:
            psi4.set_options(self.specs['psi4_options'])
        
        print("Running direct psi4 calculation")
        energy_from_psi4 , wfn_from_psi4 = psi4.opt(self.specs['method_name']+'/'+self.specs['basis_name'], molecule=self.specs['psi4_mol'], return_wfn=True)
        variables_from_psi4 = psi4.core.get_variables()
        
        print('variables from psi4')
        print(variables_from_psi4)
        
        nuclear_repulsion_energy_from_psi4=self.specs['psi4_mol'].nuclear_repulsion_energy()
        
        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()
        psi4.core.opt_clean()



        try:
            assert self.specs['psi4_options']
        except KeyError:
            my_mdt_mol.set_energy_model(mdt.models.psi4.Psi4Potential,
                                        basis=self.specs['basis_name'],
                                        memory=self.mem,
                                        theory=self.specs['method_name'])
        else:
            my_mdt_mol.set_energy_model(mdt.models.psi4.Psi4Potential,
                                        basis=self.specs['basis_name'],
                                        memory=self.mem,
                                        options=self.specs['psi4_options'],
                                        theory=self.specs['method_name'])
        
        print("Running psi4 calculation using MDT")
        nuclear_repulsion_energy_from_mdt = my_mdt_mol.minimize(requests=['opt']).__getattr__('nuclear_repulsion_energy')
                                    
        assert psi4.compare_values(nuclear_repulsion_energy_from_psi4 , nuclear_repulsion_energy_from_mdt , 6 , "Nuclear repulsion energies from optimization")

        psi4.core.clean()
        psi4.core.clean_options()
        psi4.core.clean_variables()
        psi4.core.opt_clean()


water_104_5_angle = psi4.geometry("""
    O
    H 1 0.96
    H 1 0.96 2 104.5
    """)

water_104_5_angle_freq_test = Psi4_to_Psi4_test(psi4_mol = water_104_5_angle, method_name = 'SCF', basis_name = 'cc-pVDZ')

water_104_5_angle_freq_test.freq_to_freq()


he_dimer = psi4.geometry("""
    He 0 0 0
    He 0 0 1.8
    """)

he_dimer_gradient_five_theories = [
                                   Psi4_to_Psi4_test(psi4_mol= he_dimer, method_name = 'SCF' , basis_name = 'cc-pVDZ', psi4_options = {'scf_type':'pk', 'mp2_type':'conv','reference':'rhf'}),
                                   Psi4_to_Psi4_test(psi4_mol= he_dimer, method_name = 'SCF' , basis_name = 'cc-pVTZ', psi4_options = {'scf_type':'pk', 'mp2_type':'conv','reference':'rhf'} , dertype=0) ,
                                   Psi4_to_Psi4_test(psi4_mol=he_dimer , method_name = 'SCF' , basis_name = 'cc-pV[23]Z' , psi4_options = {'scf_type':'pk', 'mp2_type':'conv','reference':'rhf'}, der_type=0 ) ,
                                   Psi4_to_Psi4_test(psi4_mol = he_dimer, method_name = 'SCF' , basis_name = 'cc-pV[23]Z' , psi4_options = {'scf_type':'pk', 'mp2_type':'conv','reference':'rhf'} ) ,
                                   Psi4_to_Psi4_test(psi4_mol = he_dimer , method_name = 'mp2' , basis_name = 'cc-pV[DT]Z' , psi4_options = {'scf_type':'pk', 'mp2_type':'conv','reference':'rhf'}) ]

for i in he_dimer_gradient_five_theories:
    i.gradient_to_gradient()

cyanic_acid = psi4.geometry("""
    H   -0.5958806528   0.9889214459   0.0000000000
    C   -0.5958806528  -0.1660941336   0.0000000000
    N    0.5535292657   0.0711607905   0.0000000000
    """)

cyanic_acid_hess_test = Psi4_to_Psi4_test(psi4_mol = cyanic_acid, method_name = 'SCF' , basis_name = 'dzp', der_type=1 )

#cyanic_acid_hess_test.hessian_to_hessian()

hydrogen = psi4.geometry("""
    H
    H 1 R
    R = 1
    """)

hydrogen_options = {'scf_type':'pk','mp2_type':'conv','g_convergence':'GAU_VERYTIGHT','e_convergence':'1.e-10'}

hydrogen_opt_tests = [Psi4_to_Psi4_test(psi4_mol = hydrogen, method_name = 'SCF', basis_name = 'cc-pVDZ' , psi4_options = hydrogen_options) ,
                      Psi4_to_Psi4_test( psi4_mol = hydrogen , method_name = 'SCF' , basis_name = 'cc-pV[DT]Z' , psi4_options = hydrogen_options ),
                      Psi4_to_Psi4_test(psi4_mol = hydrogen, method_name = 'SCF' , basis_name = 'cc-pV[DTQ]Z', psi4_options  = hydrogen_options),
                      Psi4_to_Psi4_test(psi4_mol = hydrogen , method_name = 'mp2' , basis_name = 'cc-pVDZ', psi4_options = hydrogen_options) ,
                      Psi4_to_Psi4_test(psi4_mol = hydrogen, method_name = 'mp2' , basis_name = 'cc-pV[DT]Z', psi4_options = hydrogen_options) ,
                      Psi4_to_Psi4_test(psi4_mol = hydrogen, method_name = 'mp2' , basis_name = 'cc-pV[TQ]Z' , psi4_options = hydrogen_options) ]

for i in hydrogen_opt_tests:
    i.optimize_to_optimize()


























