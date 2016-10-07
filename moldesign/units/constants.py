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

from .quantity import *

# Constants
unity = ureg.angstrom / ureg.angstrom
imi = 1.0j
pi = np.pi
sqrtpi = np.sqrt(np.pi)
sqrt2 = np.sqrt(2.0)
epsilon_0 = ureg.epsilon_0
c = ureg.speed_of_light
alpha = ureg.fine_structure_constant
hbar = ureg.hbar
boltz = boltzmann_constant = k_b = ureg.boltzmann_constant
avogadro = (1.0 * ureg.mole * ureg.avogadro_number).to(unity).magnitude

# atomic units
hartree = ureg.hartree
a0 = bohr = ureg.bohr
atomic_time = t0 = ureg.t0
electron_mass = m_e = ureg.electron_mass
electron_charge = q_e = ureg.elementary_charge

# useful units
fs = femtoseconds = ureg.fs
ps = picoseconds = ureg.ps
eV = electronvolts = ureg.eV
kcalpermol = ureg.kcalpermol
gpermol = ureg.gpermol
kjpermol = ureg.kjpermol
radians = radian = rad = ureg.rad
degrees = degree = deg = ureg.degrees
amu = da = dalton = ureg.amu
kelvin = ureg.kelvin
nm = ureg.nanometers
ang = angstrom = ureg.ang
molar = ureg.mole / ureg.liter
debye = ureg.debye

# sets default unit systems
def_length = angstrom
def_time = fs
def_vel = angstrom / fs
def_mass = amu
def_momentum = def_mass * def_vel
def_force = def_momentum / def_time
def_energy = eV