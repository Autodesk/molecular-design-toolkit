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

from moldesign import units as u
from moldesign import widgets
from . import toplevel, GeometryViewer


@toplevel
class OrbitalViewer(widgets.SelectionGroup):
    def __init__(self, mol, **kwargs):
        """
        :param mol: a molecule with A) orbitals, and B) an energy model with calculate_orbital_grid
        :param kwargs: kwargs for the viewer
        :return:
        """
        self.viewer = GeometryViewer(mol=mol, **kwargs)
        self.viewer.wfn = mol.electronic_state
        self.uipane = OrbitalUIPane(self, height=int(self.viewer.height)-50)
        hb = ipy.HBox([self.viewer, self.uipane])
        super(OrbitalViewer, self).__init__([hb])


class OrbitalUIPane(widgets.Selector, ipy.Box):
    # TODO: deal with orbitals not present in all frames of a trajectory
    # TODO: deal with orbital properties changing over a trajectory
    def __init__(self, viz, **kwargs):
        self.viz = viz
        kwargs.setdefault('width', 325)

        self.type_dropdown = ipy.Dropdown(options=self.viz.viewer.wfn.orbitals.keys())
        initial_orb = 'canonical'
        if initial_orb not in self.type_dropdown.options:
            initial_orb = self.type_dropdown.options.iterkeys().next()
        self.type_dropdown.value = initial_orb
        self.type_dropdown.observe(self.new_orb_type, 'value')

        self.orblist = ipy.Select(options={None: None},
                                  width=kwargs['width'],
                                  height=int(kwargs['height']) - 75)

        self.isoval_selector = widgets.create_value_selector(ipy.FloatSlider,
                                                          value_selects='orbital_isovalue',
                                                          min=0.0, max=0.05,
                                                          value=0.01, step=0.0001,
                                                          width=kwargs['width'],
                                                          description='Isovalue')

        self.orb_resolution = ipy.Text(description='Orbital resolution', width=75)
        self.orb_resolution.value = '40'  # this is a string to enable the 'on_submit' method
        self.change_resolution()
        self.orb_resolution.on_submit(self.change_resolution)

        children = [self.type_dropdown, self.orblist, self.isoval_selector, self.orb_resolution]
        super(OrbitalUIPane, self).__init__(children, **kwargs)
        self.new_orb_type()
        self.orblist.observe(self.new_orbital_selection, 'value')


    def new_orbital_selection(self, *args):
        self.fire_selection_event({'orbname': (self.type_dropdown.value, self.orblist.value)})

    def handle_selection_event(self, *args):
        # TODO: update the selected orbitals if something actually else triggers this
        pass

    def new_orb_type(self, *args):
        """Create list of available orbitals when user selects a new type
        """
        wfn = self.viz.viewer.wfn
        newtype = self.type_dropdown.value
        neworbs = wfn.orbitals[newtype]
        orblist = collections.OrderedDict()

        orblist[None] = None
        for i, orb in enumerate(neworbs):
            if hasattr(orb, 'unicode_name'):
                orbname = orb.unicode_name
            else:
                orbname = orb.name

            meta = ''
            if orb.energy is not None:
                meta = '{:.02fP}'.format(orb.energy.defunits())
            if orb.occupation is not None:
                if meta: meta += ', '
                meta += 'occ %.2f' % orb.occupation
            if meta:
                desc = '%d. %s   (%s)' % (i, orbname, meta)
            else:
                desc = '%d. %s' % (i, orbname)
            orblist[desc] = i
        self.orblist.value = None
        self.orblist.options = orblist

    def change_resolution(self, *args):
        viewer = self.viz.viewer
        viewer._orbital_kwargs['npts'] = int(self.orb_resolution.value)
        if viewer.current_orbital is not None:
            viewer.draw_orbital(viewer.current_orbital, render=True, **viewer._orbital_kwargs)


class VolumetricGrid(object):
    """
    Helper object for preparing gaussian CUBE files
    """
    UNITS = u.angstrom
    def __init__(self, positions, padding=2.5*u.angstrom, npoints=25):
        mins = positions.min(axis=0) - padding
        maxes = positions.max(axis=0) + padding
        self.npoints = npoints
        self.xr = (mins[0], maxes[0])
        self.yr = (mins[1], maxes[1])
        self.zr = (mins[2], maxes[2])
        self._origin = mins.value_in(self.UNITS)
        self.dx = (self.xr[1] - self.xr[0]).value_in(self.UNITS) / (float(npoints) - 1)
        self.dy = (self.yr[1] - self.yr[0]).value_in(self.UNITS) / (float(npoints) - 1)
        self.dz = (self.zr[1] - self.zr[0]).value_in(self.UNITS) / (float(npoints) - 1)
        self.fxyz = None

    def xyzlist(self):
        stride = self.npoints * 1j
        grids = np.mgrid[self.xr[0]:self.xr[1]:stride,
                self.yr[0]:self.yr[1]:stride,
                self.zr[0]:self.zr[1]:stride]
        return grids * self.UNITS

    def origin(self):
        return tuple(self._origin)

