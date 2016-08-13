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

"""
This module contains various helper functions used by MDT internally.
"""

import collections

import numpy as np
import webcolors

from moldesign import units as u


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


def get_all_atoms(*objects):
    """ Given Atoms, AtomContainers, lists of Atoms, and lists of AtomContainers,
    return a flat list of all atoms contained therein.

    A given atom is only returned once, even if it's found more than once.

    Args:
        *objects (moldesign.Atom OR moldesign.AtomContainer OR List[moldesign.Atom] OR
            List[moldesign.AtomContainer]): objects to take atoms from

    """
    from moldesign import molecules

    atoms = collections.OrderedDict()

    for obj in objects:
        if isinstance(obj, molecules.Atom):
            atoms[obj] = None
        elif hasattr(obj, 'atoms'):
            atoms.update((x,None) for x in obj.atoms)
        else:
            for item in obj:
                if isinstance(item, molecules.Atom):
                    atoms[item] = None
                elif hasattr(item, 'atoms'):
                    atoms.update((x, None) for x in item.atoms)

    return molecules.AtomList(atoms.iterkeys())


def kinetic_energy(momenta, masses):
    return 0.5 * (momenta*momenta/masses).sum()


def kinetic_temperature(ke, dof):
    from moldesign.units import k_b
    t = (2.0*ke)/(k_b*dof)
    return t.defunits()


# def get_residues(obj, **queries):
#     """
#
#     Args:
#         obj ():
#         **queries ():
#
#     Returns:
#
#     """
#     for residue in obj.residues:
#         pass
#


DEF_CATEGORICAL = 'Paired'
DEF_SEQUENTIAL = None  # should be inferno, but that's only MPL >1.5


def colormap(cats, mplmap='auto'):
    # should make it easy to choose one for:
    #  categorical data
    #  sequential (low, high important)
    #  diverging data (low, mid, high important)
    # Can deal with numerical and categorical data
    # we'll treat ints as categories for now
    global DEF_SEQUENTIAL
    from matplotlib import cm

    if hasattr(cm, 'inferno'):
        DEF_SEQUENTIAL = 'inferno'
    else:
        DEF_SEQUENTIAL = 'BrBG'

    # strip units
    units = None
    if hasattr(cats[0], 'magnitude'):
        arr = u.array(cats)
        units = arr.units
        cats = arr.magnitude

    if not isinstance(cats, np.ndarray) and not isinstance(cats[0], float):  # treat as
        # categorical
        values = np.zeros(len(cats), dtype='float')
        to_int = collections.OrderedDict()
        for i, item in enumerate(cats):
            if item not in to_int:
                to_int[item] = len(to_int)
            values[i] = to_int[item]
        if mplmap == 'auto':
            mplmap = DEF_CATEGORICAL
    else:  # it's numerical
        values = np.array(cats, dtype='float')
        if mplmap == 'auto':
            mplmap = DEF_SEQUENTIAL

    cmap = getattr(cm, mplmap)
    mx = values.max()
    mn = values.min()
    r = (values - mn) / (mx - mn)  # rescale to [0.0,1.0]
    rgb = cmap(r)
    hexcolors = [webcolors.rgb_to_hex(np.array(r[:3]) * 256) for r in rgb]
    return hexcolors