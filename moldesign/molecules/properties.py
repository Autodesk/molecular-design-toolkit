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
import numpy as np

from . import toplevel
from .. import utils


@toplevel
class MolecularProperties(utils.DotDict):
    """ Stores property values for a molecule.
    These objects will be generally created and updated by EnergyModels, not by users.
    """
    def __init__(self, mol, **properties):
        """Initialization: ``properties`` MUST include positions.

        Args:
            mol (Molecule): molecule that these properties are associated with
            **properties (dict): values of molecular properties (MUST include positions as a key)
        """
        # ADD_FEATURE: always return stored properties in the default unit systems
        positions = properties.pop('positions', mol.positions)
        super().__init__(positions=positions.copy(), **properties)

    def copy(self, mol):
        props = self.__dict__.copy()
        props['mol'] = mol
        return self.__class__(**props)

    def geometry_matches(self, mol):
        """Returns:
            bool: True if the molecule's ``position`` is the same as these properties' ``position``
        """
        return np.array_equal(self.positions, mol.positions)