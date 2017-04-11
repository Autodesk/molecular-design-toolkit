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

DEF_CATEGORICAL = 'Paired'
DEF_SEQUENTIAL = None  # should be inferno, but that's only MPL >1.5


def colormap(cats, mplmap='auto', categorical=None):
    """ Map a series of categories to hex colors, using a matplotlib colormap

    Generates both categorical and numerical colormaps.

    Args:
        cats (Iterable): list of categories or numerical values
        mplmap (str): name of matplotlib colormap object
        categorical (bool): if True, interpret this data as categorical. If False, interpret
            the data as numerical values (data must be convertible to float)

    Returns:
        List[str]: List of hexadecimal RGB color values in the in the form ``'#000102'``
    """
    # Should automatically choose the right colormaps for:
    #  categorical data
    #  sequential data (low, high important)
    #  diverging data (low, mid, high important)
    global DEF_SEQUENTIAL
    from matplotlib import cm

    if hasattr(cm, 'inferno'):
        DEF_SEQUENTIAL = 'inferno'
    else:
        DEF_SEQUENTIAL = 'BrBG'

    # strip units
    units = None  # TODO: build a color bar with units
    if hasattr(cats[0], 'magnitude'):
        arr = u.array(cats)
        units = arr.units
        cats = arr.magnitude
        is_categorical = False
    else:
        is_categorical = not isinstance(cats[0], float)

    if categorical is not None:
        is_categorical = categorical

    if is_categorical:
        values = _map_categories_to_ints(cats)
        if mplmap == 'auto':
            mplmap = DEF_CATEGORICAL
    else:
        values = np.array(map(float, cats))
        if mplmap == 'auto':
            mplmap = DEF_SEQUENTIAL

    rgb = _cmap_to_rgb(mplmap, values)
    hexcolors = [webcolors.rgb_to_hex(np.array(c)) for c in rgb]
    return hexcolors


def _map_categories_to_ints(cats):
    values = np.zeros(len(cats), dtype='float')
    to_int = collections.OrderedDict()
    for i, item in enumerate(cats):
        if item not in to_int:
            to_int[item] = len(to_int)
        values[i] = to_int[item]
    return values


def _cmap_to_rgb(mplmap, values):
    from matplotlib import cm

    cmap = getattr(cm, mplmap)
    mx = values.max()
    mn = values.min()
    cat_values = (values-mn)/(mx-mn)  # rescale values [0.0,1.0]
    rgba = cmap(cat_values)  # array of RGBA values in range [0.0, 1.0]

    # strip alpha field and rescale to [0,255] RGB integers
    rgb = [map(int, c[:3]*256.0) for c in rgba]
    return rgb


def atom_name_check(mol, force=False):
    """ Makes sure atom names are unique in each residue.

    If atoms names aren't unqiue:
      - if the names are just the names of the elements, rename them
      - else print a warning
    """
    badres = []
    for residue in mol.residues:
        names = set(atom.name for atom in residue.atoms)
        if len(names) != residue.num_atoms:
            # atom names aren't unique, check if we can change them
            for atom in residue.atoms:
                if atom.name.lower() != atom.symbol.lower():
                    badres.append(residue)
                    if not force:
                        break
            else:  # rename the atoms
                atomnums = {}
                for atom in residue.atoms:
                    atom.name = atom.symbol + str(atomnums.setdefault(atom.symbol, 0))
                    atomnums[atom.symbol] += 1

    if badres:
        print 'WARNING: residues do not have uniquely named atoms: %s' % badres
