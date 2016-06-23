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

import numpy as np
import webcolors

import moldesign as mdt
from moldesign import units as u


def get_all_atoms(*objects):
    """ Given Atoms, AtomContainers, lists of Atoms, and lists of AtomContainers,
    return a flat list of all atoms contained therein.

    A given atom is only returned once, even if it's found more than once.

    Args:
        *objects (moldesign.Atom OR moldesign.AtomContainer OR List[moldesign.Atom] OR
            List[moldesign.AtomContainer]): objects to take atoms from

    """
    atoms = collections.OrderedDict()

    for obj in objects:
        if isinstance(obj, mdt.Atom):
            atoms[obj] = None
        elif hasattr(obj, 'atoms'):
            atoms.update((x,None) for x in obj.atoms)
        else:
            for item in obj:
                if isinstance(item, mdt.Atom):
                    atoms[item] = None
                elif hasattr(item, 'atoms'):
                    atoms.update((x, None) for x in item.atoms)

    return mdt.AtomList(atoms.iterkeys())


def kinetic_energy(momenta, masses):
    return momenta.dot(momenta/(2.0*masses))


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
        arr = u.to_units_array(cats)
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