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

import traitlets
import moldesign as mdt
from moldesign import utils
from moldesign import units as u

from moldesign.utils import exports
from moldesign.uibase import ViewerToolBase, ReadoutFloatSlider

@exports
class GeometryBuilder(ViewerToolBase):
    MAXDIST = 20.0  # TODO: we need to set this dynamically
    NBR2HIGHLIGHT = '#C5AED8'
    NBR1HIGHLIGHT = '#AFC6A8'
    HIGHLIGHTOPACITY = 0.6
    POSFMT = u'{:.3f} \u212B'
    DEGFMT = u'{:.1f}\u00B0'

    def __init__(self, mol):
        super(GeometryBuilder, self).__init__(mol)

        # All numbers here are assumed angstroms and radians for now ...
        self._highlighted_bonds = []
        self._highlighted_atoms = []

        self.original_position = self.mol.positions.copy()

        self.clear_button = ipy.Button(description='Clear selection')
        self.clear_button.on_click(self.clear_selection)

        self.label_box = ipy.Checkbox(description='Label atoms', value=False)
        self.label_box.observe(self.label_atoms, 'value')

        # Viewer
        self.selection_description = ipy.HTML()
        self.subtools.children = (ipy.HBox([self.clear_button, self.label_box]),
                                  self.selection_description)
        traitlets.directional_link(
            (self.viewer, 'selected_atom_indices'),
            (self.selection_description, 'value'),
            self.get_first_atom
        )

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

        traitlets.directional_link(
            (self.viewer, 'selected_atom_indices'),
            (self.tool_holder, 'children'),
            self._get_tool_state
        )

    def get_first_atom(self, atomIndices):
        if len(atomIndices) < 1:
            return ''
        else:
            atom = self.mol.atoms[next(iter(atomIndices))]
            return u"<b>Atom</b> {atom.name} at coordinates " \
                   u"x:{p[0]:.3f}, y:{p[1]:.3f}, z:{p[2]:.3f} \u212B".format(
                       atom=atom, p=atom.position.value_in(u.angstrom))

    def set_distance(self, *args):
        bond = self.get_selected_bond(self.viewer.get_selected_bonds())
        dist_in_angstrom = self.length_slider.value
        mdt.set_distance(bond.a1, bond.a2, dist_in_angstrom*u.angstrom, adjustmol=self.adjust_button.value)
        self.viewer.set_positions()

    def set_angle(self, *args):
        bonds = self.viewer.get_selected_bonds()
        bond = self.get_selected_bond(bonds)
        bond_neighbors = self.get_bond_neighbors(bonds, bond)
        angle = self.angle_slider.value

        mdt.set_angle(bond.a1, bond.a2, bond_neighbors['a2'], angle*u.pi/180.0, adjustmol=self.adjust_button.value)
        self.viewer.set_positions()

    def set_dihedral(self, *args):
        bonds = self.viewer.get_selected_bonds()
        bond = self.get_selected_bond(bonds)
        bond_neighbors = self.get_bond_neighbors(bonds, bond)
        angle = self.dihedral_slider.value

        mdt.set_dihedral(bond_neighbors['a1'], bond.a1, bond.a2, bond_neighbors['a2'], angle*u.pi/180.0,
                                             adjustmol=self.adjust_button.value)
        self.viewer.set_positions()

    def set_atom_x(self, *args):
        pass

    def set_atom_y(self, *args):
        pass

    def set_atom_z(self, *args):
        pass

    def label_atoms(self, *args):
        self.viewer.atom_labels_shown = self.label_box.value

    # Returns the first atom indicated by atomIndices
    @staticmethod
    def get_selected_atom(atoms, atomIndices):
        return atoms[next(iter(atomIndices))]

    # Returns the first bond indicated by bondIndices
    @staticmethod
    def get_selected_bond(bonds):
        return next(iter(bonds))

    def _get_tool_state(self, selectedAtomIndices):
        # start with everything disabled
        for tool in self.atom_tools.children + self.bond_tools.children: tool.disabled = True

        if len(selectedAtomIndices) <= 0:
            return (ipy.HTML('Please click on a bond or atom'),)

        elif len(selectedAtomIndices) is 1:
            atom = self.get_selected_atom(self.mol.atoms, selectedAtomIndices)
            self.adjust_button.disabled = False

            x, y, z = atom.position.value_in(u.angstrom)
            self.x_slider.value = x
            self.x_slider.disabled = False  # for now

            self.y_slider.value = y
            self.y_slider.disabled = False  # for now

            self.z_slider.value = z
            self.z_slider.disabled = False  # for now

            for tool in self.atom_tools.children: tool.disabled = True

            return (self.atom_tools,)

        else:
            bonds = self.viewer.get_selected_bonds()
            bond = self.get_selected_bond(bonds)
            self.adjust_button.disabled = False

            # bond length
            self.length_slider.value = bond.a1.distance(bond.a2).value_in(u.angstrom)
            self.length_slider.disabled = False
            self.length_slider.description = '<b>Bond distance</b> <span style="color:{c1}">{a1.name}' \
                                             ' - {a2.name}</span>'.format(
                    a1=bond.a1, a2=bond.a2, c1=self.viewer.HIGHLIGHT_COLOR)

            a1_neighbors = set([a for a in bond.a1.bond_graph if a is not bond.a2])
            maxo = max(a1_neighbors, key=lambda x: x.mass)
            bond_neighbors = self.get_bond_neighbors(bonds, bond)

            # Bond angle
            if bond_neighbors['a2']:
                self.dihedral_slider.enable()
                self.angle_slider.enable()
                self.angle_slider.value = mdt.angle(bond.a1, bond.a2, bond_neighbors['a2']).value_in(u.degrees)
                # self.angle_slider.observe(self.set_angle, 'value')
                self.angle_slider.description = '<b>Bond angle</b> <span style="color:{c1}">{a1.name}' \
                                                ' - {a2.name}</span> ' \
                                                '- <span style="color:{c2}">{a3.name}</span>'.format(
                        a1=bond.a1, a2=bond.a2, a3=bond_neighbors['a2'],
                        c1=self.viewer.HIGHLIGHT_COLOR, c2=self.NBR2HIGHLIGHT)
            else:
                self.dihedral_slider.disable()
                self.angle_slider.disable()
                self.angle_slider.description = 'no angle associated with this bond'
                # self.angle_slider.unobserve(self.set_angle)

            # Dihedral twist
            if bond_neighbors['a2'] and bond_neighbors['a1']:
                self.dihedral_slider.value = mdt.dihedral(bond_neighbors['a1'], bond.a1, bond.a2, bond_neighbors['a2']).value_in(u.degrees)
                # self.dihedral_slider.observe(self.set_dihedral, 'value')
                self.dihedral_slider.description = '<b>Dihedral angle</b> <span style="color:{c0}">{a4.name}</span>' \
                                                   ' - <span style="color:{c1}">{a1.name}' \
                                                   ' - {a2.name}</span> ' \
                                                   '- <span style="color:{c2}">{a3.name}</span>'.format(
                        a4=bond_neighbors['a1'], a1=bond.a1, a2=bond.a2, a3=bond_neighbors['a2'],
                        c0=self.NBR1HIGHLIGHT, c1=self.viewer.HIGHLIGHT_COLOR,
                        c2=self.NBR2HIGHLIGHT)
            else:
                self.dihedral_slider.description = 'not a torsion bond'
                # self.dihedral_slider.unobserve(self.set_dihedral)
                self.dihedral_slider.disabled = True

            return [self.bond_tools]

    @staticmethod
    def get_bond_neighbors(bonds, bond):
        neighbors = { 'a1': None, 'a2': None }
        a1_neighbors = set([a for a in bond.a1.bond_graph if a is not bond.a2])
        a2_neighbors = set([a for a in bond.a2.bond_graph if a is not bond.a1])

        if a1_neighbors:
            neighbors['a1'] = max(a1_neighbors, key=lambda x: x.mass)
        if a2_neighbors:
            neighbors['a2'] = max(a2_neighbors, key=lambda x: x.mass)

        return neighbors

    def _highlight_atoms(self, atoms, color=None):
        color = utils.if_not_none(color, self.viewer.HIGHLIGHT_COLOR)
        self._highlighted_atoms += atoms
        self.viewer.add_style('vdw', atoms=atoms,
                              radius=self.viewer.ATOMRADIUS * 1.1,
                              color=color,
                              opacity=self.HIGHLIGHTOPACITY)

    def _unhighlight_atoms(self, atoms):
        self.viewer.set_style('vdw', atoms=atoms,
                              radius=self.viewer.ATOMRADIUS)

    def clear_selection(self, *args):
        self.viewer.selected_atom_indices = set()

    def reset_geometry(self, *args):
        self.clear_selection()
        self.mol.positions = self.original_position
        self.viewer.set_positions()
