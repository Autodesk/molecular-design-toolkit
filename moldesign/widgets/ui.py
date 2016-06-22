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
# TODO: catch and log event exceptions
import collections

import ipywidgets as ipy
import traitlets

import moldesign as mdt
import moldesign.geom.geometry as geo
from moldesign import utils, units as u


class Selector(object):
    """
    This is the abstract base class for something that can make a selection.
    """

    def __init__(self, *args, **kwargs):
        self.selection_group = None
        self.selection_id = None
        super(Selector, self).__init__(*args, **kwargs)

    def handle_selection_event(self, selection):
        raise NotImplementedError()

    def fire_selection_event(self, new_selection):
        self.selection_group.update_selections(self, new_selection)


class SelectionGroup(ipy.Box):
    """
    Broadcasts selections among a group of widgets.
    It doesn't do much beside rebroadcast events:
    A SelectionGroup object will call the  "handle_selection_event" methods of
    all of its children whenever its "update_selections" method is called.
    """

    def __reduce__(self):
        """These don't gat passed around,
        so it reduces to NOTHING"""
        return utils.make_none, tuple()

    def __init__(self, *args, **kwargs):
        super(SelectionGroup, self).__init__(*args, **kwargs)
        self.selection_listeners = []
        self.num_listeners = 0
        self.selection = {}
        self.viewer = None
        self.register_selection_listeners()

    def update_selections(self, event_source, selection):
        self.selection.update(selection)
        for listener in self.selection_listeners:
            listener.handle_selection_event(selection)

    def register_selection_listeners(self):
        self.num_listeners = 0
        self.selection_listeners = self.get_child_listeners(self)

    def get_child_listeners(self, element):
        if hasattr(element, 'handle_selection_event'):
            self.num_listeners += 1
            listeners = [element]
            element.selection_group = self
            element.selection_id = self.num_listeners
            if issubclass(element.__class__, mdt.viewer.GeometryViewer):
                self.viewer = element
        else:
            listeners = []

        if hasattr(element, 'children'):
            for child in element.children:
                listeners.extend(self.get_child_listeners(child))
        return listeners

    def __getattr__(self, item):
        if self.viewer is not None: return getattr(self.viewer, item)
        else: raise AttributeError(item)



class ValueSelector(Selector):
    """ This is an abstract mixin for a widget class
    """
    def __init__(self, value_selects=None, **kwargs):
        self.value_selects = value_selects
        super(ValueSelector, self).__init__(**kwargs)
        self.observe(self.value_update, 'value')
        self.__hold_fire = False

    def value_update(self, *args):
        if self.__hold_fire: return  # prevent recursive selections
        self.fire_selection_event({self.value_selects: self.value})

    def handle_selection_event(self, selection):
        self.__hold_fire = True
        try:
            if self.value_selects in selection:
                self.value = selection[self.value_selects]
        except Exception as exc:
            print 'ERROR: (ignored) %s' % exc
        self.__hold_fire = False


def create_value_selector(widget, value_selects, **kwargs):
    """
    Creates a UI element (slider, checkbox, etc.) to add to your
    selection group.
    :param widget: widget class (e.g. ipywidgets.FloatSlider)
    :param value_selects: What the value of the widget selects (e.g. time)
    :param kwargs: keyword arguments for the widget (e.g. description = "time")
    :return: the constructed widget
    """

    class SelectionWidget(ValueSelector, widget): pass

    return SelectionWidget(value_selects=value_selects, **kwargs)


class StyledTab(ipy.Tab):
    """
    Objects can inherit from this to maintain consistent styling.
    TODO: Probably better to do this with CSS?
    """

    def __init__(self, *args, **kwargs):
        kwargs.setdefault('font_size', 9)
        super(StyledTab, self).__init__(*args, **kwargs)


class AtomInspector(ipy.HTML, Selector):
    def handle_selection_event(self, selection):
        if 'atoms' not in selection: return
        atoms = selection['atoms']
        if len(atoms) == 0:
            self.value = 'No selection'
        elif len(atoms) == 1:
            atom = atoms[0]
            res = atom.residue
            chain = res.chain
            self.value = (
                "<b>Molecule</b>: %s<br>" % atom.parent.name +
                "<b>Chain</b> %s<br>" % chain.name +
                "<b>Residue</b> %s, index %d<br>" % (res.name, res.index) +
                "<b>Atom</b> %s (%s), index %d<br>" % (atom.name, atom.symbol, atom.index))
        elif len(atoms) > 1:
            atstrings = ['<b>%s</b> / res <b>%s</b> / chain <b>%s</b>' %
                         (a.name, a.residue.resname, a.chain.name)
                         for a in atoms]
            self.value = '<br>'.join(atstrings)


class ViewerToolBase(ipy.Box):
    """
    The base for most viewer-based widgets - it consists of a viewer in the top-left,
    UI controls on the right, and some additional widgets underneath the viewer
    """
    VIEWERTYPE = mdt.viewer.GeometryViewer

    def __init__(self, mol):
        self.mol = mol

        self.toolpane = ipy.Box()
        self.viewer = self.VIEWERTYPE(mol)

        self.subtools = ipy.Box()
        self.viewer_pane = ipy.VBox([self.viewer, self.subtools])
        self.main_pane = ipy.HBox([self.viewer_pane, self.toolpane])

        super(ViewerToolBase, self).__init__([self.main_pane])

    def __getattr__(self, item):
        if hasattr(self.viewer, item):
            return getattr(self.viewer, item)
        else:
            raise AttributeError(item)


class SelBase(ViewerToolBase):
    def __init__(self, mol):
        super(SelBase, self).__init__(mol)

        self._atomset = collections.OrderedDict()

        self.atom_listname = ipy.HTML('<b>Selected atoms:</b>')
        self.atom_list = ipy.SelectMultiple(options=collections.OrderedDict(),
                                            height=150)
        self.select_all_atoms_button = ipy.Button(description='Select all atoms')
        self.select_all_atoms_button.on_click(self.select_all_atoms)

        self.select_none = ipy.Button(description='Clear all selections')
        self.select_none.on_click(self.clear_selections)

        self.remove_button = ipy.Button(description='Unselect')
        self.remove_button.on_click(self.handle_remove_button_click)

    @property
    def selected_atoms(self):
        return self._atomset.keys()

    @selected_atoms.setter
    def selected_atoms(self, atoms):
        self._atomset = collections.OrderedDict((atom,None) for atom in atoms)
        self._redraw_selection_state()

    def _redraw_selection_state(self):
        self.atom_list.options = collections.OrderedDict((self.atomkey(atom), atom)
                                                         for atom in self._atomset.keys())
        self.viewer.highlight_atoms(self._atomset.keys(), render=False)
        self.viewer.render()

    def toggle_atom(self, atom):
        """Toggles atom's state in and out of the selection group"""
        if atom in self._atomset: self._atomset.pop(atom)
        else: self._atomset[atom] = None
        self._redraw_selection_state()

    def remove_atomlist_highlight(self, *args):
        self.atom_list.value = tuple()

    @staticmethod
    def atomkey(atom):
        return '%s (index %d)' % (atom.name, atom.index)

    def select_all_atoms(self, *args):
        self.selected_atoms = self.mol.atoms

    def handle_remove_button_click(self, *args):
        if self.atom_list.value:
            for atom in self.atom_list.value: self._atomset.pop(atom)
            self._redraw_selection_state()

    def clear_selections(self, *args):
        self.selected_atoms = []


class BondSelector(SelBase):
    VIEWERTYPE = mdt.viewer.BondSelectorBase

    def __init__(self, mol):
        super(BondSelector, self).__init__(mol)
        self.viewer.atom_callbacks.append(self.toggle_atom)
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


class ResidueSelector(SelBase):
    """
    Selections at the atom/residue/chain level.
    Selecting a residue selects all of its atoms.
    Selecting all atoms of a residue is equivalent to selecting the residue.
    A residue is not selected if only some of its atoms are selected.
    """
    def __init__(self, mol):
        super(ResidueSelector, self).__init__(mol)
        self.viewer.add_click_callback(self.atom_click)

        self._residue_selection = collections.OrderedDict()
        self._residueset = collections.OrderedDict()
        self.selection_type = ipy.Dropdown(description='Clicks select:',value='Residue',
                                           options=('Atom', 'Residue', 'Chain'))

        self.residue_listname = ipy.HTML('<b>Selected residues:</b>')
        self.residue_list = ipy.SelectMultiple(options=collections.OrderedDict(),
                                               height=150)
        self.residue_list.observe(self.remove_atomlist_highlight, 'value')
        self.atom_list.observe(self.remove_reslist_highlight, 'value')

        self.subtools.children = [ipy.HBox([self.select_all_atoms_button, self.select_none])]
        self.toolpane.children = [self.selection_type,
                                  self.atom_listname,
                                  self.atom_list,
                                  self.residue_listname,
                                  self.residue_list,
                                  self.remove_button]

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

    def atom_click(self, atom):
        if self.selection_type.value == 'Atom':
            self.toggle_atom(atom)
        elif self.selection_type.value == 'Residue':
            self.toggle_residue(atom.residue, clickatom=atom)
        elif self.selection_type.value == 'Chain':
            self.toggle_chain(atom.chain, atom)
        else:
            raise ValueError('Unknown selecton_type %s' % self.selection_type.value)

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

class ReadoutFloatSlider(ipy.Box):
    description = traitlets.Unicode()
    value = traitlets.Float()

    def __init__(self, format=None, *args, **kwargs):
        description = kwargs.pop('description', 'FloatSlider')
        min = kwargs.setdefault('min', 0.0)
        max = kwargs.setdefault('max', 10.0)
        self.formatstring = format
        self.header = ipy.HTML()
        self.readout = ipy.Text(width=100)
        self.readout.on_submit(self.parse_value)

        kwargs.setdefault('readout', False)
        self.slider = ipy.FloatSlider(*args, **kwargs)
        self.minlabel = ipy.HTML(u'<font size=1.5>{}</font>'.format(self.formatstring.format(min)))
        self.maxlabel = ipy.HTML(u'<font size=1.5>{}</font>'.format(self.formatstring.format(max)))
        self.sliderbox = ipy.HBox([self.minlabel, self.slider, self.maxlabel])
        traitlets.link((self, 'description'), (self.header, 'value'))
        traitlets.link((self, 'value'), (self.slider, 'value'))
        self.description = description
        self.update_readout()
        super(ReadoutFloatSlider, self).__init__([self.header,
                                                  self.readout,
                                                  self.sliderbox])

    @traitlets.observe('value')
    def update_readout(self, *args):
        self.readout.value = self.formatstring.format(self.value)

    def parse_value(self, *args):
        try:
            f = float(self.readout.value)
        except ValueError:
            s = self.readout.value
            match = utils.GETFLOAT.search(s)
            if match is None:
                self.readout.value = self.formatstring.format(self.slider.value)
                print "Couldn't parse string %s" % s
                return
            else:
                f = float(s[match.start():match.end()])
        self.slider.value = f


class GeometryBuilder(ViewerToolBase):
    VIEWERTYPE = mdt.viewer.BondSelectorBase

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
        geo.set_distance(sel.a1, sel.a2, dist_in_angstrom * u.angstrom, adjustmol=self.adjust_button.value)
        self.viewer.set_positions()

    def set_angle(self, *args):
        sel = self._selection
        assert sel.type == 'bond'
        angle = self.angle_slider.value
        geo.set_angle(sel.a1, sel.a2, sel.nbr_a2, angle * u.pi / 180.0, adjustmol=self.adjust_button.value)
        self.viewer.set_positions()

    def set_dihedral(self, *args):
        sel = self._selection
        assert sel.type == 'bond'
        angle = self.dihedral_slider.value
        geo.set_dihedral(sel.nbr_a1, sel.a1, sel.a2, sel.nbr_a2, angle * u.pi / 180.0,
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
                self.angle_slider.value = geo.angle(sel.a1, sel.a2, sel.nbr_a2).value_in(u.degrees)
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
                self.dihedral_slider.value = geo.dihedral(sel.nbr_a1, sel.a1, sel.a2, sel.nbr_a2).value_in(u.degrees)
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
