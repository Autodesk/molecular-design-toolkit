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
import ipywidgets as ipy

import moldesign as mdt
from moldesign import utils
from moldesign.viewer import BondClicker
from moldesign import units as u

from moldesign.uibase import ViewerToolBase, ReadoutFloatSlider

def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class GeometryBuilder(ViewerToolBase):
    VIEWERTYPE = BondClicker

    MAXDIST = 20.0  # TODO: we need to set this dynamically
    NBR2HIGHLIGHT = '#C5AED8'
    NBR1HIGHLIGHT = '#AFC6A8'
    HIGHLIGHTOPACITY = 0.6
    POSFMT = u'{:.3f} \u212B'
    DEGFMT = u'{:.1f}\u00B0'

    def __init__(self, mol):
        super(GeometryBuilder, self).__init__(mol)

        # All numbers here are assumed angstroms and radians for now ...
        self._selection = utils.DotDict(blank=True, type=None)
        self._highlighted_bonds = []
        self._highlighted_atoms = []

        self.original_position = self.mol.positions.copy()

        self.clear_button = ipy.Button(description='Clear selection')
        self.clear_button.on_click(self.clear_selection)

        self.label_box = ipy.Checkbox(description='Label atoms', value=False)
        self.label_box.observe(self.label_atoms, 'value')

        # Viewer
        self.viewer.atom_callbacks.append(self.atom_click)
        self.viewer.bond_callbacks.append(self.bond_click)

        self.selection_description = ipy.HTML()

        self.subtools.children = (ipy.HBox([self.clear_button, self.label_box]),
                                  self.selection_description)

        # Atom manipulation tools
        self.x_slider = ReadoutFloatSlider(min=-self.MAXDIST, max=self.MAXDIST,
                                           description='<b>x</b>', format=self.POSFMT)
        self.x_slider.observe(self.set_atom_x, 'value')
        self.y_slider = ReadoutFloatSlider(min=-self.MAXDIST, max=self.MAXDIST,
                                           description='<b>y</b>', format=self.POSFMT)
        self.y_slider.observe(self.set_atom_y, 'value')

        self.z_slider = ReadoutFloatSlider(min=-self.MAXDIST, max=self.MAXDIST,
                                           description='<b>z</b>', format=self.POSFMT)
        self.z_slider.observe(self.set_atom_z, 'value')

        # Bond manipulation tools
        self.adjust_button = ipy.Checkbox(description='Adjust entire molecule', align='end', value=True)

        self.length_slider = ReadoutFloatSlider(min=0.1, max=self.MAXDIST, format=self.POSFMT)
        self.length_slider.observe(self.set_distance, 'value')
        self.angle_slider = ReadoutFloatSlider(min=1.0, max=179.0, step=2.0, format=self.DEGFMT)
        self.angle_slider.observe(self.set_angle, 'value')
        self.dihedral_slider = ReadoutFloatSlider(min=-90.0, max=360.0, step=4.0, format=self.DEGFMT)
        self.dihedral_slider.observe(self.set_dihedral, 'value')

        self.bond_tools = ipy.VBox((self.adjust_button,
                                    self.length_slider,
                                    self.angle_slider,
                                    self.dihedral_slider))

        self.atom_tools = ipy.VBox((self.adjust_button,
                                    self.x_slider,
                                    self.y_slider,
                                    self.z_slider))

        self.reset_button = ipy.Button(description='Reset geometry')
        self.reset_button.on_click(self.reset_geometry)

        self.tool_holder = ipy.VBox()
        self.toolpane.children = (self.tool_holder,
                                  self.reset_button)

    def set_distance(self, *args):
        sel = self._selection
        assert sel.type == 'bond'
        dist_in_angstrom = self.length_slider.value
        mdt.set_distance(sel.a1, sel.a2, dist_in_angstrom*u.angstrom, adjustmol=self.adjust_button.value)
        self.viewer.set_positions()

    def set_angle(self, *args):
        sel = self._selection
        assert sel.type == 'bond'
        angle = self.angle_slider.value
        mdt.set_angle(sel.a1, sel.a2, sel.nbr_a2, angle*u.pi/180.0, adjustmol=self.adjust_button.value)
        self.viewer.set_positions()

    def set_dihedral(self, *args):
        sel = self._selection
        assert sel.type == 'bond'
        angle = self.dihedral_slider.value
        mdt.set_dihedral(sel.nbr_a1, sel.a1, sel.a2, sel.nbr_a2, angle*u.pi/180.0,
                                             adjustmol=self.adjust_button.value)
        self.viewer.set_positions()

    def set_atom_x(self, *args):
        pass

    def set_atom_y(self, *args):
        pass

    def set_atom_z(self, *args):
        pass

    def label_atoms(self, *args):
        if self.label_box.value:
            self.viewer.label_atoms()
        else:
            self.viewer.remove_all_labels()

    def atom_click(self, atom):
        sel = self._selection
        if sel.blank:  # select this atom
            sel.blank = False
            sel.type = 'atom'
            sel.atom = atom

        elif sel.type == 'atom':  # We've selected 2 atoms - i.e. a bond
            if atom is sel.atom:  # clicked twice -> deselect the thing
                return self.clear_selection()

            elif atom in sel.atom.bond_graph:  # select the bond
                return self.bond_click(mdt.Bond(sel.atom, atom))  # turn this into a bond selection
            else:  # select a new atom
                self.clear_selection(render=False)
                sel = self._selection
                sel.blank = False
                sel.type = 'atom'
                sel.atom = atom

        elif sel.type == 'bond':
            if atom in sel.a1_neighbors:  # change the neighboring selection
                sel.nbr_a1 = atom
            elif atom in sel.a2_neighbors:
                sel.nbr_a2 = atom
            else:  # select a new atom
                self.clear_selection(render=False)
                return self.atom_click(atom)

        self._redraw_selection()

    def _set_tool_state(self):
        # start with everything disabled
        for tool in self.atom_tools.children + self.bond_tools.children: tool.disabled = True

        if self._selection.blank:
            self.tool_holder.children = (ipy.HTML('Please click on a bond or atom'),)

        elif self._selection.type == 'atom':
            self.adjust_button.disabled = False

            self.tool_holder.children = (self.atom_tools,)
            x, y, z = self._selection.atom.position.value_in(u.angstrom)
            self.x_slider.value = x
            self.x_slider.disabled = False  # for now

            self.y_slider.value = y
            self.y_slider.disabled = False  # for now

            self.z_slider.value = z
            self.z_slider.disabled = False  # for now

            for tool in self.atom_tools.children: tool.disabled = True

        elif self._selection.type == 'bond':
            sel = self._selection
            self.adjust_button.disabled = False

            # bond length
            self.length_slider.value = sel.a1.distance(sel.a2).value_in(u.angstrom)
            self.length_slider.disabled = False
            self.length_slider.description = '<b>Bond distance</b> <span style="color:{c1}">{a1.name}' \
                                             ' - {a2.name}</span>'.format(
                    a1=sel.a1, a2=sel.a2, c1=self.viewer.HIGHLIGHT_COLOR)

            # Bond angle
            if sel.nbr_a2:
                self.angle_slider.value = mdt.angle(sel.a1, sel.a2, sel.nbr_a2).value_in(u.degrees)
                # self.angle_slider.observe(self.set_angle, 'value')
                self.angle_slider.disabled = False
                self.angle_slider.description = '<b>Bond angle</b> <span style="color:{c1}">{a1.name}' \
                                                ' - {a2.name}</span> ' \
                                                '- <span style="color:{c2}">{a3.name}</span>'.format(
                        a1=sel.a1, a2=sel.a2, a3=sel.nbr_a2,
                        c1=self.viewer.HIGHLIGHT_COLOR, c2=self.NBR2HIGHLIGHT)
            else:
                self.angle_slider.description = 'no angle associated with this bond'
                # self.angle_slider.unobserve(self.set_angle)
                self.angle_slider.disabled = True

            # Dihedral twist
            if sel.nbr_a2 and sel.nbr_a1:
                self.dihedral_slider.value = mdt.dihedral(sel.nbr_a1, sel.a1, sel.a2, sel.nbr_a2).value_in(u.degrees)
                # self.dihedral_slider.observe(self.set_dihedral, 'value')
                self.dihedral_slider.disabled = False
                self.dihedral_slider.description = '<b>Dihedral angle</b> <span style="color:{c0}">{a4.name}</span>' \
                                                   ' - <span style="color:{c1}">{a1.name}' \
                                                   ' - {a2.name}</span> ' \
                                                   '- <span style="color:{c2}">{a3.name}</span>'.format(
                        a4=sel.nbr_a1, a1=sel.a1, a2=sel.a2, a3=sel.nbr_a2,
                        c0=self.NBR1HIGHLIGHT, c1=self.viewer.HIGHLIGHT_COLOR,
                        c2=self.NBR2HIGHLIGHT)
            else:
                self.dihedral_slider.description = 'not a torsion bond'
                # self.dihedral_slider.unobserve(self.set_dihedral)
                self.dihedral_slider.disabled = True

            self.tool_holder.children = [self.bond_tools]

        else:
            raise ValueError('Unknown selection type %s' % self._selection.type)

    def bond_click(self, bond):
        sel = self._selection
        if sel.type == 'bond':  # check if this bond is already selected
            a1, a2 = bond.a1, bond.a2
            if (a1 is sel.a1 and a2 is sel.a2) or (a1 is sel.a2 and a2 is sel.a1):
                return self.clear_selection()

        self.clear_selection(render=False)
        sel = self._selection
        sel.blank = False
        sel.type = 'bond'
        sel.bond = bond
        sel.a1 = bond.a1
        sel.a2 = bond.a2
        sel.a1_neighbors = set([a for a in bond.a1.bond_graph if a is not bond.a2])
        sel.a2_neighbors = set([a for a in bond.a2.bond_graph if a is not bond.a1])
        sel.nbr_a1 = sel.nbr_a2 = None
        if sel.a1_neighbors:
            sel.nbr_a1 = max(sel.a1_neighbors, key=lambda x: x.mass)
        if sel.a2_neighbors:
            sel.nbr_a2 = max(sel.a2_neighbors, key=lambda x: x.mass)
        self._redraw_selection()

    def _highlight_atoms(self, atoms, color=None, render=True):
        color = utils.if_not_none(color, self.viewer.HIGHLIGHT_COLOR)
        self._highlighted_atoms += atoms
        self.viewer.add_style('vdw', atoms=atoms,
                              radius=self.viewer.ATOMRADIUS * 1.1,
                              color=color,
                              opacity=self.HIGHLIGHTOPACITY,
                              render=render)

    def _unhighlight_atoms(self, atoms, render=True):
        self.viewer.set_style('vdw', atoms=atoms,
                              radius=self.viewer.ATOMRADIUS,
                              render=render)

    def _redraw_selection(self):
        # unhighlight any previous selections
        if self._highlighted_atoms:
            self._unhighlight_atoms(self._highlighted_atoms, render=False)
        self._highlighted_atoms = []
        for bond in self._highlighted_bonds:
            self.viewer.unset_bond_color(bond, render=False)
        self._highlighted_bonds = []

        # Set the selection view
        sel = self._selection
        if sel.type == 'atom':
            self._highlight_atoms([sel.atom], render=False)
            self.selection_description.value = \
                u"<b>Atom</b> {atom.name} at coordinates " \
                u"x:{p[0]:.3f}, y:{p[1]:.3f}, z:{p[2]:.3f} \u212B".format(
                        atom=sel.atom, p=sel.atom.position.value_in(u.angstrom))

        elif sel.type == 'bond':
            self.selection_description.value = "<b>Bond:</b> %s - %s" % (sel.a1.name, sel.a2.name)
            self._highlighted_bonds = [sel.bond]
            self.viewer.set_bond_color(self.viewer.HIGHLIGHT_COLOR, sel.bond, render=False)
            self._highlight_atoms([sel.a1, sel.a2], render=False)

            if sel.nbr_a1 is not None:
                nmdtond = mdt.Bond(sel.a1, sel.nbr_a1)
                self._highlight_atoms([sel.nbr_a1], color=self.NBR1HIGHLIGHT, render=False)
                self.viewer.set_bond_color(self.NBR1HIGHLIGHT, nmdtond, render=False)
                self._highlighted_bonds.append(nmdtond)
            if sel.nbr_a2 is not None:
                nmdtond = mdt.Bond(sel.a2, sel.nbr_a2)
                self._highlight_atoms([sel.nbr_a2], color=self.NBR2HIGHLIGHT, render=False)
                self.viewer.set_bond_color(self.NBR2HIGHLIGHT, nmdtond, render=False)
                self._highlighted_bonds.append(nmdtond)

        elif sel.type is not None:
            raise ValueError('Unknown selection type %s' % self._selection.type)

        self.viewer.render()
        self._set_tool_state()

    def clear_selection(self, render=True, *args):
        self._selection = utils.DotDict(blank=True, type=None)
        self.selection_description.value = ""
        if render: self._redraw_selection()

    def reset_geometry(self, *args):
        self.clear_selection(render=False)
        self.mol.positions = self.original_position
        self.viewer.set_positions()
        self._redraw_selection()