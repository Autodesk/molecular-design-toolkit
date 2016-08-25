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

from moldesign import viewer
from moldesign.uibase.components import SelBase


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class BondSelector(SelBase):
    VIEWERTYPE = viewer.BondClicker

    def __init__(self, mol):
        super(BondSelector, self).__init__(mol)
        self.viewer.bond_callbacks.append(self.toggle_bond)

        self._bondset = collections.OrderedDict()
        self._drawn_bond_state = set()

        self.bond_listname = ipy.HTML('<b>Selected bonds:</b>')
        self.bond_list = ipy.SelectMultiple(options=collections.OrderedDict(),
                                            height=150)
        self.bond_list.observe(self.remove_atomlist_highlight, 'value')
        self.atom_list.observe(self.remove_bondlist_highlight, 'value')

        self.select_all_bonds_button = ipy.Button(description='Select all bonds')
        self.select_all_bonds_button.on_click(self.select_all_bonds)

        self.subtools.children = [ipy.HBox([self.select_all_atoms_button,
                                            self.select_all_bonds_button,
                                            self.select_none])]
        self.toolpane.children = (self.atom_listname,
                                  self.atom_list,
                                  self.bond_listname,
                                  self.bond_list,
                                  self.remove_button)

    def select_all_bonds(self, *args):
        self.selected_bonds = list(self.mol.bonds)

    @property
    def selected_bonds(self):
        return self._bondset.keys()

    @selected_bonds.setter
    def selected_bonds(self, newbonds):
        self._bondset = collections.OrderedDict((b,None) for b in newbonds)
        self._redraw_selection_state()

    def _redraw_selection_state(self):
        currentset = set(self._bondset)

        to_turn_on = currentset.difference(self._drawn_bond_state)
        to_turn_off = self._drawn_bond_state.difference(currentset)
        for bond in to_turn_off: self.viewer.unset_bond_color(bond, render=False)
        for b in to_turn_on: self.viewer.set_bond_color(self.viewer.HIGHLIGHT_COLOR, b, render=False)

        self.bond_list.options = collections.OrderedDict((self.bondkey(bond), bond) for bond in self._bondset)
        super(BondSelector, self)._redraw_selection_state()
        self._drawn_bond_state = currentset
        self.remove_atomlist_highlight()

    def remove_bondlist_highlight(self, *args):
        self.bond_list.value = tuple()

    @staticmethod
    def bondkey(bond):
        return bond.name

    def toggle_bond(self, bond):
        if bond in self._bondset: self._bondset.pop(bond)  # unselect the bond
        else: self._bondset[bond] = None  # select the bond
        self._redraw_selection_state()

    def handle_remove_button_click(self, *args):
        if self.bond_list.value:
            for bond in self.bond_list.value: self._bondset.pop(bond)
            self._redraw_selection_state()
        super(BondSelector, self).handle_remove_button_click(*args)

    def select_all_bonds(self, *args):
        self.selected_bonds = list(self.mol.bonds)

    def clear_selections(self, *args):
        self.selected_bonds = []
        super(BondSelector, self).clear_selections(*args)


@exports
class ResidueSelector(SelBase):
    """
    Selections at the atom/residue/chain level.
    Selecting a residue selects all of its atoms.
    Selecting all atoms of a residue is equivalent to selecting the residue.
    A residue is not selected if only some of its atoms are selected.
    """
    def __init__(self, mol):
        super(ResidueSelector, self).__init__(mol)

        self._residue_selection = collections.OrderedDict()
        self._residueset = collections.OrderedDict()
        self.selection_type = ipy.Dropdown(description='Clicks select:',value=self.viewer.selection_type,
                                           options=('Atom', 'Residue', 'Chain'))

        traitlets.link((self.selection_type, 'value'), (self.viewer, 'selection_type'))

        self.residue_listname = ipy.HTML('<b>Selected residues:</b>')
        self.residue_list = ipy.SelectMultiple(options=list(), height=150)
        traitlets.directional_link((self.viewer, 'selected_atoms'), (self.residue_list, 'options'), self._atoms_to_residues)
        self.residue_list.observe(self.remove_atomlist_highlight, 'value')
        self.atom_list.observe(self.remove_reslist_highlight, 'value')

        self.subtools.children = [ipy.HBox([self.select_all_atoms_button, self.select_none])]
        self.toolpane.children = [self.selection_type,
                                  self.atom_listname,
                                  self.atom_list,
                                  self.residue_listname,
                                  self.residue_list,
                                  self.remove_button]

    # Returns a list of all residues indicated in the set of atoms
    def _atoms_to_residues(self, selected_atoms):
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
        for atomIndex in selected_atoms:
            atom = self.mol.atoms[atomIndex]
            if 'selected_count' in residues[atom.residue.index]:
                residues[atom.residue.index]['selected_count'] += 1
            else:
                residues[atom.residue.index]['selected_count'] = 1

            if residues[atom.residue.index]['selected_count'] >= residues[atom.residue.index]['total']:
                selected_residues.add(atom.residue.name)

        return list(selected_residues)

    def _redraw_selection_state(self):
        # this is slow and crappy ...
        super(ResidueSelector, self)._redraw_selection_state()

        # Update the residue list
        def pop_residue(r):
            resopts.pop(self.reskey(r))
            self._residueset.pop(r)

        resopts = self.residue_list.options.copy()
        atomcounts = collections.Counter()
        for atom in self._atomset: atomcounts[atom.residue] += 1
        for res in atomcounts:
            if res.num_atoms == atomcounts[res]:  # i.e., this residue IS fully selected
                if res not in self._residueset:
                    resopts[self.reskey(res)] = res
                    self._residueset[res] = None
            else:  # i.e., this residue should NOT be selected
                if res in self._residueset: pop_residue(res)

        for res in self._residueset:
            if res not in atomcounts:
                pop_residue(res)

        self.residue_list.options = resopts

    @property
    def selected_residues(self):
        return self._residueset.keys()

    @selected_residues.setter
    def selected_residues(self, residues):
        newres = set(residues)

        for res in newres.symmetric_difference(self._residueset):
            self.toggle_residue(res, render=False)

        self._residueset = newres
        self._redraw_selection_state()

    def toggle_residue(self, residue, clickatom=None, render=True):
        if clickatom is not None:
            deselect = (clickatom in self._atomset)
        else:
            deselect = (residue in self._residueset)

        if deselect:
            for atom in residue.atoms: self._atomset.pop(atom, None)
        else:
            for atom in residue.atoms: self._atomset[atom] = None

        if render: self._redraw_selection_state()

    def remove_reslist_highlight(self, *args):
        self.atom_list.value = tuple()

    @staticmethod
    def atomkey(atom):
        return '%s (index %d)' % (atom.name, atom.index)

    @staticmethod
    def reskey(residue):
        return '{res.name} in chain "{res.chain.name}"'.format(res=residue)

    def handle_remove_button_click(self, *args):
        if self.residue_list.value:
            for res in self.residue_list.value:
                self.toggle_residue(res, render=False)
            self._redraw_selection_state()
        else:
            super(ResidueSelector, self).handle_remove_button_click(*args)
