from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library
standard_library.install_aliases()

# Copyright 2017 Autodesk Inc.
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
import numpy as np

import moldesign as mdt
import moldesign.molecules.atomcollections
from moldesign import external
from moldesign.mathutils import sub_angles, apply_4x4_transform

from . import toplevel, angle, dihedral
from .coords import _infer_dihedral


@toplevel
def set_distance(a1, a2, newlength, adjustmol=True):
    """ Set the distance between two atoms. They will be adjusted along the vector separating them.
    If the two atoms are A) bonded, B) not part of the same ring system, and C) ``adjustmol`` is
    True, then the entire molecule's positions will be modified as well

    Args:
        a1,a2 (mdt.Atom): atoms to adjust
        newlength (u.Scalar[length]): new length to set
        adjustmol (bool): Adjust all atoms on either side of this bond?
    """
    #TODO: lots of room for optimization here
    if adjustmol:
        assert a1.molecule is not None
        assert a1.molecule == a2.molecule
    vec = a1.position - a2.position
    dist = np.sqrt(vec.dot(vec))
    direction = (vec / dist)
    delta = newlength - dist
    if np.abs(delta) < 1.0e-5 * delta.get_units(): return
    if not adjustmol:
        a1.position += direction * delta / 2.0
        a2.position -= direction * delta / 2.0
    else:
        mol = a1.molecule
        indices, sign = _get_fragment_indices(mol, a1, a2)
        mol.positions[indices] += delta*direction*sign


@toplevel
def set_angle(a1, a2, a3, theta, adjustmol=True):
    """ Set the angle between bonds a1-a2 and a3-a2. The atoms will be adjusted along the
    gradient of the angle.
    If ``adjustmol`` is True and the topology is unambiguous, then the entire molecule's positions
    will be modified as well

    Args:
        a1,a2,a3 (mdt.Atom): atoms to adjust
        theta (u.Scalar[angle]): new angle to set
        adjustmol (bool): Adjust all atoms on either side of this bond?
    """
    # TODO: deal with co-linear a1, a2, a3 - the rotation axis is ill-defined in this case \
    #      (require an axis to be specified)
    # TODO: weakly cache the rotation axis so that users can set angle to 0 or 180 without losing the axis
    current = angle(a1, a2, a3)
    rotation = sub_angles(current, theta)
    if abs(rotation) < 1.0e-6: return

    axis = np.cross(a1.position - a2.position, a3.position - a2.position) # do vecs need to be normalized?
    if not adjustmol:
        rotmat_l = external.transformations.rotation_matrix(rotation / 2.0, axis, a2.position)
        rotmat_r = external.transformations.rotation_matrix(-rotation / 2.0, axis, a2.position)

        a1.position = apply_4x4_transform(rotmat_l, a1.position)
        a3.position = apply_4x4_transform(rotmat_r, a3.position)

    else:
        mol = a2.molecule
        indices, sign = _get_fragment_indices(mol, a1, a2)
        rotmat = external.transformations.rotation_matrix(rotation, axis*sign, a2.position)
        mol.positions[indices] = apply_4x4_transform(rotmat, mol.positions[indices])


@toplevel
def set_dihedral(a1, a2=None, a3=None, a4=None, theta=None, adjustmol=True):
    """ Set the twist angle of atoms a1 and a4 around the central bond a2-a3. The atoms will be
    adjusted along the gradient of the angle.

    Can be called as ``set_dihedral(a1, a2, a3, a4, theta, adjustmol=True)``
              OR     ``set_dihedral(a2, a2, theta, adjustmol=True)``
              OR     ``set_dihedral(bond, theta, adjustmol=True)``

    If ``adjustmol`` is True and the topology is unambiguous, then the entire molecule's positions
    will be modified as well

    Args:
        a1 (mdt.Bond): central bond in dihedral
        a1,a2 (mdt.Atom): atoms around central bond in dihedral
        a3, a4 (mdt.Atom):
        theta (u.Scalar[angle]): new angle to set
        adjustmol (bool): Adjust all atoms on either side of this bond?
    """
    # TODO: deal with co-linear a1/a4, a2, a3 - the angle is ill-defined \
    #      (should just an arbitrary axis normal to the central bond)
    if a4 is None:
        if isinstance(a1, mdt.Bond):
            if theta is None:
                theta = a2
            a1, a2 = a1.a1, a1.a2
        if a3 is not None and theta is None:
            theta, a3 = a3, theta
        elif a3 is not None or a4 is not None or theta is None:
            raise ValueError('Invalid number of arguments for set_dihedral')
        a1, a2, a3, a4 = _infer_dihedral(a1, a2)

    current = dihedral(a1, a2, a3, a4)
    rotation = sub_angles(theta, current)
    if abs(rotation) < 1.0e-6: return

    axis = a2.position - a3.position
    if not adjustmol:
        rotmat_l = external.transformations.rotation_matrix((-rotation / 2.0), axis, a3.position)
        rotmat_r = external.transformations.rotation_matrix((rotation / 2.0), axis, a3.position)

        a1.position = apply_4x4_transform(rotmat_l, a1.position)
        a4.position = apply_4x4_transform(rotmat_r, a4.position)

    else:
        mol = a2.molecule
        indices, sign = _get_fragment_indices(mol, a3, a2)
        rotmat = external.transformations.rotation_matrix(rotation, axis*sign, a3.position)
        mol.positions[indices] = apply_4x4_transform(rotmat, mol.positions[indices])


def _get_fragment(mol, a1, a2):
    """
    Given a pair of atoms a1 and a2, return two fragments, one composed of all atoms on a1's
    side of the molecule, the other composed of all atoms on a2's side of the molecule.

    This won't work if a1 and a2 are in a cycle.
    """
    # DFS for a1's side of the bond. To prevent visiting a2, we mark it as visited at the start
    visited = set([a2])
    def dfs_dive(atom):
        visited.add(atom)
        for nbr in atom.bond_graph:
            if nbr is a2 and atom is not a1:
                raise ValueError("a1 and a2 are in a cyclic moiety")
            if nbr not in visited:
                dfs_dive(nbr)
    dfs_dive(a1)
    visited.remove(a2)
    result = moldesign.molecules.atomcollections.AtomList(visited)
    return result


def _get_fragment_indices(mol, a1, a2):
    key = (mol, a1, a2)
    if key in _get_fragment_indices.cache:
        return _get_fragment_indices.cache[key]

    # Try to get the smaller fragment ... a bit hacky right now
    try: frag1 = _get_fragment(mol, a1, a2)
    except ValueError:
        frag1 = None

    try: frag2 = _get_fragment(mol, a2, a1)
    except ValueError:
        if frag1 is None: raise
        else:
            frag = frag1
            sign = 1.0
    else:
        if frag1 is None or len(frag1) > len(frag2):
            frag = frag2
            sign = -1.0
        else:
            frag = frag1
            sign = 1.0

    indices = [atom.index for atom in frag]
    result = np.array(indices)
    _get_fragment_indices.cache[key] = (result, sign)
    return result, sign

_get_fragment_indices.cache = {}
