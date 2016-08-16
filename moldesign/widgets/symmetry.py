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
import numpy as np

import moldesign as mdt
from moldesign import units as u


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class Symmetrizer(ipy.Box):
    def __init__(self, mol):
        self._current_shapes = []
        self.mol = mol
        self.tolerance = 0.3 * u.angstrom
        self.original_coords = mol.positions.copy()

        self.showing = ipy.HTML()
        self.viewer = mol.draw3d()
        """:type viewer: moldesign.viewer.GeometryViewer"""

        self.description = ipy.HTML()
        self.symm_selector = ipy.Select()
        self.symm_selector.observe(self.show_symmetry, names='value')

        self.apply_button = ipy.Button(description='Symmetrize')
        self.apply_button.on_click(self.apply_selected_symmetry)

        self.reset_button = ipy.Button(description='Reset')
        self.reset_button.on_click(self.reset_coords)

        self.apply_all_button = ipy.Button(description='Apply all', padding=10)
        self.apply_all_button.on_click(self.set_highest_symmetry)

        self.tolerance_descrip = ipy.HTML(u'<small>tolerance/\u212B</small>',)
        self.tolerance_chooser = ipy.BoundedFloatText(value=self.tolerance.value_in(u.angstrom),
                                                      min=0.0)

        self.recalculate_button = ipy.Button(description='Recalculate')
        self.recalculate_button.on_click(self.coords_changed)

        self.symm_pane = ipy.VBox([self.description,
                                   self.symm_selector,
                                   ipy.HBox([self.apply_button, self.reset_button]),
                                   self.apply_all_button,
                                   ipy.HBox([self.tolerance_chooser, self.recalculate_button]),
                                   self.tolerance_descrip],
                                  width=325)

        self.symmetry = None
        self.coords_changed()

        self.hbox = ipy.HBox([ipy.VBox([self.viewer, self.showing]), self.symm_pane])
        super(Symmetrizer, self).__init__([self.hbox])

    def reset_coords(self, *args):
        self.mol.positions = self.original_coords
        self.viewer.append_frame(positions=self.original_coords)
        self.coords_changed()

    def coords_changed(self, *args):
        with self.symm_selector.hold_trait_notifications():
            self.symm_selector.options = {}
        self.description.value = 'Finding symmetries ...'
        self.tolerance = self.tolerance_chooser.value * u.angstrom
        self.symmetry = mdt.geom.get_symmetry(self.mol, tolerance=self.tolerance)
        options = collections.OrderedDict()
        for elem in self.symmetry.elems:
            if elem.max_diff.magnitude != 0.0:
                key = '{0}) {1} (error={2.magnitude:.4f} {2.units})'.format(elem.idx, elem.symbol, elem.max_diff)
            else:
                key = '{0}) {1} (exact)'.format(elem.idx, elem.symbol, elem.max_diff)
            options[key] = elem
        with self.symm_selector.hold_trait_notifications():
            self.symm_selector.options = options
        descrip = 'Highest symmetry group: <b>%s</b><br>' % self.symmetry.symbol
        if self.symmetry.rms.magnitude == 0.0:
            descrip += '(Exact)'
        else:
            descrip += 'RMS Error = {:.03P}'.format(self.symmetry.rms)
        self.description.value = descrip
        self.viewer.append_frame(positions=self.symmetry.orientation)
        self.viewer.center()

    def apply_selected_symmetry(self, *args):
        idx = self.symm_selector.value.idx
        elem = self.symmetry.elems[idx]
        newcoords = self.symmetry.get_symmetrized_coords(elem)
        self.mol.atoms.position = newcoords
        if not np.allclose(newcoords, self.symmetry.orientation, atol=1.0e-10):
            self.viewer.append_frame(positions=newcoords)
            self.coords_changed()

    def show_symmetry(self, *args):
        self.showing.value = ''
        if self._current_shapes:
            for s in self._current_shapes: self.viewer.remove(s, render=False)
            self._current_shapes = []
        if self.symm_selector.value is None:
            self.viewer.render()
            return

        elem = self.symm_selector.value
        symbol = elem.symbol

        self.showing.value = '%s visualization not implemented' % symbol

        if symbol == 'C1':
            self.viewer.render()
            self.showing.value = 'Identity operation'
            return

        elif symbol == 'Ci':
            inversion = self.viewer.draw_sphere(np.zeros(3) * u.angstrom,
                                                color='0x4AB4C4',
                                                radius=0.5 * u.angstrom,
                                                opacity=0.85, render=False)
            self._current_shapes.append(inversion)
            self.showing.value = 'Inversion center'

        elif symbol == 'Cs' or (symbol[0] == 'S' and symbol[1].isdigit()):
            axis = elem.get_axis()
            rad = 2.5 * max(self.symmetry.orientation.max(), 3.0 * u.angstrom)
            plane = self.viewer.draw_circle(np.zeros(3),
                                            axis,
                                            radius=rad,
                                            opacity=0.6,
                                            color='0xAB00FE', render=False)
            self._current_shapes.append(plane)
            self.showing.value = 'Mirror plane (normal = %s)' % axis

        if symbol[0] in 'SC' and symbol[1].isdigit():
            axis = elem.get_axis()
            nrot = int(symbol[1])
            projections = self.symmetry.orientation.dot(axis)
            top = axis * max(3.25 * projections.max(), 3.0*u.angstrom)
            bottom = axis * min(2.5 * projections.min(), -2.5*u.angstrom)
            arrow = self.viewer.draw_arrow(start=bottom, end=top,
                                           color='0x00FE03', render=False, opacity=0.8)
            self._current_shapes.append(arrow)
            if symbol[0] == 'S':
                self.showing.value = '%d-fold improper rotation axis (%s)' % (nrot, axis)
            else:
                self.showing.value = '%d-fold rotation axis (%s)' % (nrot, axis)

        self.viewer.render()


    def set_highest_symmetry(self, *args):
        raise NotImplementedError()