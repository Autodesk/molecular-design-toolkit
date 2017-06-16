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
from ..compute import runsremotely
from ..interfaces import psi4_interface
from .base import QMBase

# On June 15, JSOB added " import numpy as np
import numpy as np

# On June 13, JSOB added " ,'mp2':'mp2' "

THEORIES = {'rhf':'scf', 'mp2':'mp2'}

BASES = {}

# NATIVE_OPTIONS added by JSOB
NATIVE_OPTIONS = {'mdt_frozen_core':'freeze_core','you_name_this_one':'how_about_not'}

class Psi4Potential(QMBase):
    def prep(self):
        return True

    @runsremotely(enable=psi4_interface.force_remote)
    def calculate(self, requests):
        import psi4
        
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

        geom = psi4_interface.mdt_to_psi4(self.mol)
        self.name = psi4_interface.mol_name(self.mol)
                        #molecule=geom added on June 13 by sherrills
        
        #output=self._get_runmethod(requests)(self._gettheorystring(), molecule=geom,return_wfn=True)
        #props={'potential_energy':output[0] * u.hartree, 'psi4_wavefunction':output[1]}
        
        
        #props = {'potential_energy': (self._get_runmethod(requests)(self._gettheorystring(), molecule=geom,return_wfn=True))[0]} * u.hartree
        #            psi4.energy                  ('mp2')
        props=self._get_properties(requests, geom)

        psi4.core.clean()
        psi4.core.clean_options()
        return props

    @staticmethod
    def _get_runmethod(requests):
        import psi4
        # TODO: actually deal with requests
        return psi4.energy

    def _gettheorystring(self):
        return THEORIES.get(self.params.theory, self.params.theory)

    def _getoptionsstring(self):
        for mdt_opt in self.params.options:
            if mdt_opt in NATIVE_OPTIONS:
                self.params.options[NATIVE_OPTIONS[mdt_opt]]=self.params.options[mdt_opt]
                del self.params.options[mdt_opt]
        return self.params.options

    def _get_properties(self,requests,geom):
        import psi4
        props={}
#WRITE FOR THIS FUNCTION OR OTHER FUNCTIONS APPROPRIATELY IN RELATION TO RETURN STATEMENT
        psi4_output = self._get_runmethod(requests)(self._gettheorystring(), molecule=geom,return_wfn=True)
        
        props['potential_energy']=psi4_output[0]
        
# i.e. psi4_output[0] is the scf energy psi4_output[1] is the wavefunction

# Find oeprop() elements that is "One-electron quantities" as recognized in psi4
        psi4_vars_title=self.name+'_psi4_get_variables'
        psi4.oeprop(psi4_output[1], 'DIPOLE',title=psi4_vars_title)
        psi4_variables=psi4.core.get_variables()
                #tmp_array=[ psi4_variables[psi4_vars_title+' DIPOLE X'], psi4_variables[psi4_vars_title+' DIPOLE Y'],psi4_variables[psi4_vars_title+' DIPOLE Z'])
        #placeholder= psi4_vars_title+' DIPOLE'
        #props['dipole'] = [ X, Y, Z for keys in psi4_vars_title  ]
        #props['dipole'] = [ psi4_variables[psi4_vars_title+' DIPOLE X'],
        #                    psi4_variables[psi4_vars_title+' DIPOLE Y'],
        #                    psi4_variables[psi4_vars_title+' DIPOLE Z'] ]
        
        props['dipole'] = [ psi4_variables[x] for x in list(psi4_variables.keys()) if x[:-2] == psi4_vars_title+' DIPOLE']



        return props
########################################################
# DO NOT END _get_properties ABOVE OR BELOW THIS LINE
########################################################

