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

from moldesign import units as u

__all__ = 'ATOMIC_MASSES ATOMIC_NUMBERS ELEMENTS SYMBOLS'.split()

ATOMIC_NUMBERS = {'Ac': 89, 'Ag': 47, 'Al': 13, 'Am': 95, 'Ar': 18, 'As': 33, 'At': 85, 'Au': 79,
                  'B': 5, 'Ba': 56, 'Be': 4, 'Bh': 107, 'Bi': 83, 'Bk': 97, 'Br': 35, 'C': 6,
                  'Ca': 20, 'Cd': 48, 'Ce': 58, 'Cf': 98, 'Cl': 17, 'Cm': 96, 'Cn': 112, 'Co': 27,
                  'Cr': 24, 'Cs': 55, 'Cu': 29, 'Db': 105, 'Ds': 110, 'Dy': 66, 'Er': 68, 'Es': 99,
                  'Eu': 63, 'F': 9, 'Fe': 26, 'Fm': 100, 'Fr': 87, 'Ga': 31, 'Gd': 64, 'Ge': 32,
                  'H': 1, 'He': 2, 'Hf': 72, 'Hg': 80, 'Ho': 67, 'Hs': 108, 'I': 53, 'In': 49, 'Ir': 77,
                  'K': 19, 'Kr': 36, 'La': 57, 'Li': 3, 'Lr': 103, 'Lu': 71, 'Md': 101, 'Mg': 12,
                  'Mn': 25, 'Mo': 42, 'Mt': 109, 'N': 7, 'Na': 11, 'Nb': 41, 'Nd': 60, 'Ne': 10,
                  'Ni': 28, 'No': 102, 'Np': 93, 'O': 8, 'Os': 76, 'P': 15, 'Pa': 91, 'Pb': 82,
                  'Pd': 46, 'Pm': 61, 'Po': 84, 'Pr': 59, 'Pt': 78, 'Pu': 94, 'Ra': 88, 'Rb': 37,
                  'Re': 75, 'Rf': 104, 'Rg': 111, 'Rh': 45, 'Rn': 86, 'Ru': 44, 'S': 16, 'Sb': 51,
                  'Sc': 21, 'Se': 34, 'Sg': 106, 'Si': 14, 'Sm': 62, 'Sn': 50, 'Sr': 38, 'Ta': 73,
                  'Tb': 65, 'Tc': 43, 'Te': 52, 'Th': 90, 'Ti': 22, 'Tl': 81, 'Tm': 69, 'U': 92,
                  'Uuh': 116, 'Uuo': 118, 'Uup': 115, 'Uuq': 114, 'Uus': 117, 'Uut': 113, 'V': 23,
                  'W': 74, 'Xe': 54, 'Y': 39, 'Yb': 70, 'Zn': 30, 'Zr': 40}
ELEMENTS = {val: key for key, val in ATOMIC_NUMBERS.items()}
SYMBOLS = ELEMENTS


# Isotopic masses for the most abundant species of each element
# from https://www.ncsu.edu/chemistry/msf/pdf/IsotopicMass_NaturalAbundance.pdf
ATOMIC_MASSES = {i: m*u.amu for i, m in zip(range(1, 55), (
    1.007825, 4.002603, 7.016004, 9.012182, 11.009305, 12.0, 14.003074, 15.994915, 18.998403,
    19.99244, 22.98977,
    23.985042, 26.981538, 27.976927, 30.973762, 31.972071, 34.968853, 39.962383, 38.963707,
    39.962591, 44.95591,
    47.947947, 50.943964, 51.940512, 54.93805, 55.934942, 58.9332, 57.935348, 62.929601, 63.929147,
    68.925581,
    73.921178, 74.921596, 79.916522, 78.918338, 83.911507, 84.911789, 87.905614, 88.905848,
    89.904704, 92.906378,
    97.905408, 97.907216, 101.90435, 102.905504, 107.903894, 106.905093, 113.903358, 114.903878,
    119.902197, 120.903818,
    129.906223, 126.904468, 131.904154))}

for atnum, mass in list(ATOMIC_MASSES.items()):
    ATOMIC_MASSES[ELEMENTS[atnum]] = mass  # index by atnum and symbol

ATOMIC_MASSES[-1] = -1.0*u.amu
