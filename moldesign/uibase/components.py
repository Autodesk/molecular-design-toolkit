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

from moldesign import utils, viewer
from moldesign import units as u

from .selector import Selector


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
            lines = ["<b>Molecule</b>: %s<br>" % atom.molecule.name]
            if atom.chain.name is not None:
                lines.append("<b>Chain</b> %s<br>" % chain.name)
            if atom.residue.type != 'placeholder':
                lines.append("<b>Residue</b> %s, index %d<br>" % (res.name, res.index))
            lines.append("<b>Atom</b> %s (%s), index %d<br>" % (atom.name, atom.symbol, atom.index))
            self.value = '\n'.join(lines)

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
    VIEWERTYPE = viewer.GeometryViewer

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


class ReadOnlyRepr(ipy.Box):
    """ When a value is assigned, displays its __repr__ instead
    """
    def __init__(self, *args, **kwargs):
        super(ReadOnlyRepr, self).__init__(*args, **kwargs)
        self.textbox = ipy.Text()
        self.textbox.disabled = True
        self.children = [self.textbox]

    @property
    def value(self):
        return self.textbox.value

    @value.setter
    def value(self, v):
        self.textbox.value = repr(v)


class UnitText(ipy.Box):
    """Widget for a user to input a quantity with physical units.

    Args:
        value (u.Scalar): initial value for the widget
        units (u.MdtUnit): if provided, require that input has the same dimensionality as these
            units

    Note:
        A very rough approximation of the FloatText, IntText, etc. widgets - this does NOT have a
        counterpart in JavaScript. But it will behave similarly when setting units and values
        from python.

    """
    INVALID = u'\u274C'
    VALID = u"\u2705"

    def __init__(self, value=None, units=None, **kwargs):
        super(UnitText, self).__init__(layout=ipy.Layout(display='flex', flex_flow='row wrap'),
                                       **kwargs)
        self.textbox = ipy.Text()
        self.textbox.observe(self._validate, 'value')
        self._error_msg = None

        if units is not None:
            self.dimensionality = u.get_units(units).dimensionality
        else:
            self.dimensionality = None

        self._validated_value = None
        self.validated = ipy.HTML(self.INVALID)
        self.children = [self.textbox, self.validated]
        self._is_valid = False
        if value is not None:
            self.value = value

    def _validate(self, change):
        self._validated_value = None
        self._error_msg = False

        # Check that we can parse this
        try:
            val = u.ureg(change['new'])
        except:  # parsing failed, pint just raises generic "Exception"
            self._error_msg = "Failed to parse '%s'" % self.textbox.value
        else:  # Check dimensionality
            valdim = u.get_units(val).dimensionality
            if self.dimensionality is not None and valdim != self.dimensionality:
                self._error_msg = "Requires dimensionality %s" % self.dimensionality

        if not self._error_msg:
            self.validated.value = '<span title="Valid quantity">%s</span>' % self.VALID
            self._validated_value = val
        else:
            self.validated.value = '<span title="%s">%s</span>' % (self._error_msg, self.INVALID)

    @property
    def value(self):
        if self._validated_value is None:
            raise ValueError(self._error_msg)
        else:
            return self._validated_value

    @value.setter
    def value(self, v):
        self.textbox.value = str(v)
