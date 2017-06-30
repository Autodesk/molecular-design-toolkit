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

import os
import yaml
from moldesign import units as u

from . import PACKAGEPATH

__all__ = 'ATOMIC_MASSES ATOMIC_NUMBERS ELEMENTS SYMBOLS'.split()

with open(os.path.join(PACKAGEPATH, '_static_data', 'nist_atomic.yml'), 'r') as ymlfile:
    isotopes = yaml.load(ymlfile)

ATOMIC_NUMBERS = {record[0]['symbol']: atnum for atnum, record in isotopes.items()}
ELEMENTS = {val: key for key, val in ATOMIC_NUMBERS.items()}
SYMBOLS = ELEMENTS


# Isotopic masses for the most abundant species of each element
# from https://www.ncsu.edu/chemistry/msf/pdf/IsotopicMass_NaturalAbundance.pdf
ATOMIC_MASSES = {atnum: records[0]['mass']*u.amu for atnum, records in isotopes.items()}


for atnum, mass in list(ATOMIC_MASSES.items()):
    ATOMIC_MASSES[ELEMENTS[atnum]] = mass  # index by atnum and symbol

ATOMIC_MASSES[-1] = -1.0*u.amu
