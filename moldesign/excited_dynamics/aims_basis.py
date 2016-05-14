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
from moldesign.excited_dynamics import wavepackets
from moldesign.base import propagators
from moldesign.units import *
from base.atom import Atom


class AIMS_Atom(Atom):
    def __init__(self,*args,**kwargs):
        super(AIMS_Atom,self).__init__(*args,**kwargs)
        self.width = 0.0*def_length


class AIMS_BasisFn(propagators.VelocityVerletMD, wavepackets.GaussianWavepacket):
    """
    Holds a time-dependent nuclear/electronic product basis function.
    In practical terms, this is a nuclear configuration with
    specified gaussian widths, an associated electronic wavefunction.
    """
    def __init__(self,*args,**kwargs):
        super(AIMS_BasisFn,self).__init__(*args,**kwargs)
        self.state_id=-1
        self.semiclassical_phase = 0.0 * radians
        self.amplitude = 0.0j





class AIMS_Wavefn(object):
    """
    Holds a time-dependent nuclear/electronic product basis function.
    In practical terms, this is a nuclear configuration with
    specified gaussian widths, an associated electronic wavefunction.
    """
    def __init__(self):
        pass
