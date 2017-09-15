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
from . import toplevel
from .. import utils, units


@toplevel
class MolecularProperties(utils.DotDict):
    """ Stores property values for a molecule.
    These objects will be generally created and updated by EnergyModels, not by users.

    Args:
        mol (Molecule): molecule that these properties are associated with
        **properties (dict): values of molecular properties (MUST include positions as a key)
    """
    def __init__(self, mol, **properties):
        self.mol = mol
        positions = properties.pop('positions', mol.positions)
        super().__init__(positions=positions.copy(), **properties)

    def copy_to(self, mol):
        """
        Args:
            mol (moldesign.Molecule): molecule to copy these properties to

        Returns:
            MolecularProperties: copied instance, associated with the new molecule
        """
        newprops = self.__class__(mol, positions=self.positions)
        for name, val in self.items():
            if name == 'positions':
                continue
            elif hasattr(val, 'copy_to'):
                newprops[name] = val.copy_to(mol)
            elif hasattr(val, 'copy'):
                newprops[name] = val.copy()
            else:
                newprops[name] = val
        return newprops

    def geometry_matches(self, mol):
        """ Returns true if a molecule's positions are the same as in these properties.

        Allows for some very slight numerical noise due to units conversion and associated issues.

        Returns:
            bool: True if the molecule's ``position`` is the same as these properties' ``position``
        """
        return units.arrays_almost_equal(self.positions, mol.positions)

    def __getitem__(self, item):
        return units.default.convert_if_possible(super().__getitem__(item))

    def __getattr__(self, item):
        return units.default.convert_if_possible(super().__getattr__(item))


class AtomicProperties(utils.NewUserDict):
    """
    Stores calculated atomic properties

    Internally, references atoms by their indices rather than by object, which lets us
    copy molecular properties without weird reference issues

    Args:
        mapping (Mapping[moldesign.Atom, Any]): map of atoms to properties
    """
    type = 'atomic'  # backwards-compatible signal that this is a map from atoms to property values

    def __init__(self, mapping=None):
        super().__init__()
        if mapping is not None:
            self.update(mapping)

    def copy(self):
        new = self.__class__()
        super(self.__class__, new).update(self)
        return new

    def update(self, mapping):
        for k,v in mapping.items():
            self[k] = v

    # TODO: check that passed atoms actually belong to the correct molecule
    def __setitem__(self, atom, value):
        super().__setitem__(self._getkey(atom), value)

    def __getitem__(self, atom):
        return super().__getitem__(self._getkey(atom))

    @staticmethod
    def _getkey(atom):
        if isinstance(atom, int):
            k = atom
        else:
            k = atom.index
        return k

    def __contains__(self, atom):
        return super().__contains__(self._getkey(atom))
