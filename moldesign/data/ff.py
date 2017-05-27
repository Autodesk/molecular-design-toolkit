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

from . import PACKAGEPATH
from ..utils import exports_names

exports_names('AMBER_DEFAULT', 'AMBER_LEAPRC', 'AMBER_SYSTEM')

AMBER_DEFAULT = ('GLYCAM_06j-1', 'tip3p', 'gaff2', 'lipid14', 'OL15', 'OL3', 'ff14SB')

with open(os.path.join(PACKAGEPATH, '_static_data', 'amber_ffs.yml'), 'r') as _ambfile:
    ambff = yaml.load(_ambfile)

AMBER_LEAPRC = {}
AMBER_SYSTEM = {}
for system, ffs in list(ambff.items()):
    for ff, leaprc in list(ffs.items()):
        AMBER_LEAPRC[ff] = leaprc
        AMBER_SYSTEM[ff] = system
