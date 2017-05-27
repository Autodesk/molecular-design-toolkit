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
#from singledispatch import singledispatch
import numpy as np

from moldesign import units as u
#from moldesign.molecules import Bond
from moldesign.mathutils import normalized, safe_arccos

from . import toplevel


# NEWFEATURE: preliminary profiling indicates that, for interactive work, the UNITS LIBRARY is
#       actually the biggest bottleneck

# NEWFEATURE: pyramidalization aka out-of-plane bending

@toplevel
#@singledispatch
def distance(a1, a2):
    """ Return distance between two atoms

    Args:
        a1,a2 (mdt.Atom): the two atoms

    Returns:
        u.Scalar[length]: the distance
    """
    diffvec = a1.position - a2.position
    dim = len(diffvec.shape)
    if dim == 1:
        return np.sqrt(diffvec.dot(diffvec))
    else:
        assert dim == 2
        return np.sqrt((diffvec*diffvec).sum(axis=1))


# @distance.register(Bond)
# def distance(bond):
#     return distance(bond.a1, bond.a2)


@toplevel
#@singledispatch
def angle(a1, a2, a3):
    """ The angle between bonds a2-a1 and a2-a3

    Args:
        a1,a2,a3 (mdt.Atom): the atoms describing the angle

    Returns:
        u.Scalar[length]: the distance
    """
    r21 = (a1.position - a2.position).defunits_value()  # remove units immediately to improve speed
    r23 = (a3.position - a2.position).defunits_value()
    e12 = normalized(r21)
    e23 = normalized(r23)

    dim = len(e12.shape)
    if dim == 2:
        costheta = (e12*e23).sum(axis=1)
    else:
        assert dim == 1
        costheta = np.dot(e12, e23)
    theta = safe_arccos(costheta)
    return theta * u.radians

#
# @angle.register(Bond)
# def angle(b1, b2):
#     a1, a2, a3 = _join_bonds(b1, b2)
#     return angle(a1, a2, a3)


def _join_bonds(b1, b2):
    """ Return a1, a2, a3, where a2 is the atom shared by both bonds
    """
    try:
        a3 = b2.partner(b1.a2)
    except ValueError:
        return b1.a2, b1.a1, b2.partner(b1.a1)
    else:
        return b1.a1, b1.a2, a3


@toplevel
#@singledispatch
def dihedral(a1, a2=None, a3=None, a4=None):
    """Twist angle of bonds a1-a2 and a4-a3 around around the central bond a2-a3

    Can be called as ``dihedral(a1, a2, a3, a4)``
              OR     ``dihedral(a2, a2)``
              OR     ``dihedral(bond)``

    Args:
        a1 (mdt.Bond): the central bond in the dihedral. OR
        a1,a2 (mdt.Atom): the atoms describing the dihedral
        a3,a4 (mdt.Atom): (optional) if not passed, ``a1`` and ``a2`` will be treated as the
            central atoms in this bond, and a3 and a4 will be inferred.

    Returns:
        (units.Scalar[angle]): angle -  [0, 2 pi) radians
    """
    if a3 is a4 is None:  # infer the first and last atoms
        a1, a2, a3, a4 = _infer_dihedral(a1, a2)

    r21 = (a1.position - a2.position).defunits_value()  # remove units immediately to improve speed
    r34 = (a4.position - a3.position).defunits_value()

    dim = len(r21.shape)
    center_bond = (a2.position - a3.position).defunits_value()
    plane_normal = normalized(center_bond)

    if dim == 1:
        r21_proj = r21 - plane_normal * r21.dot(plane_normal)
        r34_proj = r34 - plane_normal * r34.dot(plane_normal)
    else:
        assert dim == 2
        r21_proj = r21 - plane_normal * (r21*plane_normal).sum(axis=1)[:, None]
        r34_proj = r34 - plane_normal * (r34*plane_normal).sum(axis=1)[:, None]

    va = normalized(r21_proj)
    vb = normalized(r34_proj)

    if dim == 1:
        costheta = np.dot(va, vb)
        if np.allclose(costheta, 1.0):
            return 0.0 * u.radians
        elif np.allclose(costheta, -1.0):
            return u.pi * u.radians
    else:
        costheta = (va*vb).sum(axis=1)

    abstheta = safe_arccos(costheta)
    cross = np.cross(va, vb)

    if dim == 1:
        theta = abstheta * np.sign(np.dot(cross, plane_normal))
    else:
        theta = abstheta * np.sign((cross*plane_normal).sum(axis=1))
    return (theta * u.radians) % (2.0 * u.pi * u.radians)


def _infer_dihedral(a2, a3=None):
    """ Given two atoms defining the central bond in a dihedral, pick the first and last atoms
    in a heuristic way (see :meth:`_pick_atom`) to get a unique-ish definition.
    """
    if a3 is None:  # assume bond-like
        bond = a2
        a2, a3 = bond.a1, bond.a2
    a1 = _pick_atom(a2, a3)
    a4 = _pick_atom(a3, a2)
    return a1, a2, a3, a4


def _pick_atom(atom, nbr):
    """ Pick an atom bonded to ``atom`` that:
      A) is not nbr
      B) has the largest atomic number
      C) has the lowest index

    This gives us a unique definition for dihedrals when only passing 2 atoms
    """
    newat = None

    if hasattr(atom, 'traj'):
        istraj = True
        traj = atom.traj
        atom = atom.real_atom
        nbr = nbr.real_atom
    else:
        istraj = False

    for bond in atom.bonds:
        pt = bond.partner(atom)
        if pt == nbr:
            continue
        elif newat is None:
            newat = pt
        elif newat.atnum < pt.atnum:
            newat = pt
        elif newat.atnum == pt.atnum and newat.index > pt.index:
            newat = pt

    if newat is None:
        raise ValueError('%s is not part of a dihedral' % atom)

    if istraj:
        newat = traj.atoms[newat.index]

    return newat




#
# @dihedral.register(Bond)
# def dihedral(b1, b2=None, b3=None):
#     if b2 is b3 is None:
#         return dihedral(b1.a1, b1.a2)
#
#     else:
#         a1, a2, a3 = _join_bonds(b1, b2)
#         a22, a33, a4 = _join_bonds(b2, b3)
#
#         assert a22 == a2
#         assert a33 == a3
#
#         return dihedral(a1, a2, a3, a4)
#
#
