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
import collections

import ipywidgets as ipy
import traitlets

from moldesign import viewer, utils
from moldesign.uibase.components import SelBase


@utils.exports
class BondSelector(SelBase):
    def __init__(self, mol):
        super(BondSelector, self).__init__(mol)

        self._bondset = collections.OrderedDict()
        self._drawn_bond_state = set()

        self.bond_listname = ipy.HTML('<b>Selected bonds:</b>')
        self.bond_list = ipy.SelectMultiple(options=list(),
                                            height=150)

        traitlets.directional_link(
            (self.viewer, 'selected_atom_indices'),
            (self.bond_list, 'options'),
            self._atoms_to_bonds
        )

        self.atom_list.observe(self.remove_bondlist_highlight, 'value')

        self.subtools.children = [ipy.HBox([self.select_all_atoms_button,
                                            self.select_none])]
        self.toolpane.children = (self.atom_listname,
                                  self.atom_list,
                                  self.bond_listname,
                                  self.bond_list)

    def _atoms_to_bonds(self, atomIndices):
        return list(self.viewer.get_selected_bonds(atomIndices))

    def _redraw_selection_state(self):
        currentset = set(self._bondset)

        to_turn_on = currentset.difference(self._drawn_bond_state)
        to_turn_off = self._drawn_bond_state.difference(currentset)

        self.bond_list.options = collections.OrderedDict((self.bondkey(bond), bond) for bond in self._bondset)
        super(BondSelector, self)._redraw_selection_state()
        self._drawn_bond_state = currentset
        self.remove_atomlist_highlight()

    def remove_bondlist_highlight(self, *args):
        self.bond_list.value = tuple()

    @staticmethod
    def bondkey(bond):
        return bond.name

    def clear_selections(self, *args):
        super(BondSelector, self).clear_selections(*args)


@utils.exports
class ResidueSelector(SelBase):
    """
    Selections at the atom/residue/chain level.
    Selecting a residue selects all of its atoms.
    Selecting all atoms of a residue is equivalent to selecting the residue.
    A residue is not selected if only some of its atoms are selected.
    """
    def __init__(self, mol):
        super(ResidueSelector, self).__init__(mol)

        self.selection_type = ipy.Dropdown(description='Clicks select:',value=self.viewer.selection_type,
                                           options=('Atom', 'Residue', 'Chain'))

        traitlets.link((self.selection_type, 'value'), (self.viewer, 'selection_type'))

        self.residue_listname = ipy.HTML('<b>Selected residues:</b>')
        self.residue_list = ipy.SelectMultiple(options=list(), height=150)
        traitlets.directional_link((self.viewer, 'selected_atom_indices'), (self.residue_list, 'options'), self._atoms_to_residues)
        self.residue_list.observe(self.remove_atomlist_highlight, 'value')
        self.atom_list.observe(self.remove_reslist_highlight, 'value')

        self.subtools.children = [ipy.HBox([self.select_all_atoms_button, self.select_none])]
        self.toolpane.children = [self.selection_type,
                                  self.atom_listname,
                                  self.atom_list,
                                  self.residue_listname,
                                  self.residue_list]

    # Returns a list of all residues indicated in the set of atoms
    def _atoms_to_residues(self, selected_atom_indices):
        # Get all residues and their total number of atoms
        residues = dict()
        for atom in self.mol.atoms:
            if atom.residue.index in residues:
                residues[atom.residue.index]['total'] += 1
            else:
                residues[atom.residue.index] = {
                    'total': 1,
                }

        # Get the total number of selected atoms per residue and compare to total
        selected_residues = set()
        for atomIndex in selected_atom_indices:
            atom = self.mol.atoms[atomIndex]
            if 'selected_count' in residues[atom.residue.index]:
                residues[atom.residue.index]['selected_count'] += 1
            else:
                residues[atom.residue.index]['selected_count'] = 1

            if residues[atom.residue.index]['selected_count'] >= residues[atom.residue.index]['total']:
                selected_residues.add(atom.residue)

        return list(selected_residues)

    @property
    def selected_residues(self):
        return self._atoms_to_residues(self.viewer.selected_atom_indices)

    @selected_residues.setter
    def selected_residues(self, residues):
        self.viewer.select_residues(residues)

    def toggle_residue(self, residue):
        self.viewer.toggle_residues([residue])

    def remove_reslist_highlight(self, *args):
        self.atom_list.value = tuple()

    @staticmethod
    def atomkey(atom):
        return '%s (index %d)' % (atom.name, atom.index)

    @staticmethod
    def reskey(residue):
        return '{res.name} in chain "{res.chain.name}"'.format(res=residue)
