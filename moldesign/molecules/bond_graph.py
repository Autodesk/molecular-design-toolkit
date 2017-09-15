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
from future.utils import PY2
import collections

from .. import utils

if not PY2:  # we'll skip the type hints in Python 2
    from typing import Mapping
    import moldesign as mdt  # needed for the late-binding type hints


class AtomBondDict(utils.NewUserDict):
    """ Holds the bonds for a given atom

    This object does the legwork of making sure that the bond graph remains symmetric - it makes
    sure any modifications to this data are reflected on the symmetric side of the graph
    """
    if not PY2:
        __metaclass__ = Mapping['mdt.Atom', int]

    def __init__(self, parent, atom):
        super().__init__()
        self.parent = parent
        self.atom = atom

    def __setitem__(self, otheratom, order):
        self._setitem_super(otheratom, order)
        if self.parent is not None:
            self.parent[otheratom]._setitem_super(self.atom, order)

    def _setitem_super(self, atom, order):
        super().__setitem__(atom, order)

    def __delitem__(self, otheratom):
        self._delitem_super(otheratom)
        if self.parent is not None:
            self.parent[otheratom]._delitem_super(self.atom)

    def _delitem_super(self, otheratom):
        super().__delitem__(otheratom)


class BondGraph(utils.NewUserDict):
    """ The bond graph of a molecule. This object manages the molecular bonding topology.

    This object can generally be treated as a dict, with two important exceptions:
      1. it maintains symmetry, so that setting ``bond_graph[a][b] = o`` will automatically set
         ``bond_graph[b][a] = o``
      2. It has exactly one entry for every atom in the molecule; it does not expose any methods
         for adding or deleting these entries.
    """
    if not PY2:
        __metaclass__ = Mapping['mdt.Atom', AtomBondDict]

    def __init__(self, molecule):
        """ Initializes the object, and raises up the bond state of all atoms in the molecule
        """
        super().__init__()
        self.molecule = molecule
        self._add_atoms(molecule.atoms)

    def __str__(self):
        'Bond graph for %s' % self.molecule

    def __repr__(self):
        return '<%s>' % self.molecule

    def _add_atoms(self, atoms):
        """ Move bonds defined in the atom to this object.

        If the atom is already present in the graph, do nothing.

        This method raises state - it centralizes bonding information about the molecule in
        one place. However, bonding information about atoms OUTSIDE of the molecule is left
        untouched.

        Args:
            atoms (List[moldesign.Atom]): atom to initialize
        """
        add_atoms = [atom for atom in atoms if atom not in self]
        for atom in add_atoms:
            super().__setitem__(atom, AtomBondDict(self, atom))
        for atom in add_atoms:
            bg = atom._bond_graph
            for nbr, order in list(bg.items()):
                if nbr in self:
                    self[atom][nbr] = order
                    del bg[nbr]

    def __setitem__(self, key, args):
        raise NotImplementedError('Bond graphs cannot be manipulated using this method')

    def __delitem__(self, key):
        raise NotImplementedError('Bond graphs cannot be manipulated using this method')



