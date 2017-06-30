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
from future.moves.collections import UserDict

import moldesign as mdt  # needed for the late-binding type hints
from typing import Mapping


class _DNull(object):
    pass


class AtomBondDict(UserDict, Mapping['mdt.Atom', int]):
    """ Holds the bonds for a given atom
    """
    def __init__(self, parent, atom):
        super().__init__()
        self.parent = parent
        self.atom = atom

    def __setitem__(self, otheratom, order):
        super().__setitem__(otheratom, order)
        if self.parent is not None:
            otherdict = self.parent[otheratom]
            super(otherdict.__class__, otherdict).__setitem__(self.atom, order)

    def pop(self, otheratom, default=_DNull):
        if otheratom in self:
            retval = self[otheratom]
            super().__delitem__(otheratom)
            if self.parent is not None:
                otherdict = self.parent[otheratom]
                super(otherdict.__class__, otherdict).__delitem__(self.atom)
            return retval

        elif default is not _DNull:
            return default

        else:
            raise KeyError("'%s' is not a key in %s" % self)

    def __delitem__(self, otheratom):
        self.pop(otheratom)


class BondGraph(UserDict, Mapping['mdt.Atom', AtomBondDict]):
    """ The bond graph for this molecule.

    This object can generally be treated as a dict, with two important exceptions:
      1. it maintains symmetry, so that setting ``bond_graph[a][b] = o`` will automatically set
         ``bond_graph[b][a] = o``
      2. It has exactly one entry for every atom in the molecule; it does not expose any methods
         for adding or deleting these entries.
    """
    def __init__(self, molecule):
        """ Initializes the object, and raises up the bond state of all atoms in the molecule
        """
        super().__init__()
        self.molecule = molecule
        self._add_atoms(molecule.atoms)

    def _add_atoms(self, atoms):
        """ Takes control of this atom's bonds

        If the atom is already present in the graph, nothing happens

        Args:
            atoms (List[moldesign.Atom]): atom to initialize
        """
        add_atoms = [atom for atom in atoms if atom not in self]
        for atom in add_atoms:
            super().__setitem__(atom, AtomBondDict(self, atom))
        for atom in add_atoms:
            bg = atom._bond_graph
            for nbr, order in bg.items():
                if nbr in self:
                    self[atom][nbr] = order
            atom._bond_graph = None

    def __setitem__(self, key, args):
        raise NotImplementedError('Bond graphs cannot be manipulated using this method')

    def pop(self, *args):
        raise NotImplementedError('Bond graphs cannot be manipulated using this method')

    def __delitem__(self, key):
        raise NotImplementedError('Bond graphs cannot be manipulated using this method')



