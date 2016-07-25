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
import numpy as np

from moldesign import units as u
from moldesign.mathutils import normalized, safe_arccos

from . import toplevel


# NEWFEATURE: preliminary profiling indicates that, for interactive work, the UNITS LIBRARY is
#       actually the biggest bottleneck

# NEWFEATURE: pyramidalization aka out-of-plane bending

@toplevel
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
        return np.sqrt((diffvec*diffvec).sum(axis=0))


@toplevel
def angle(a1, a2, a3):
    """ The angle between bonds a2-a1 and a2-a3
    Args:
        a1,a2,a3 (mdt.Atom): the atoms describing the vector

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


@toplevel
def dihedral(a1, a2, a3, a4):
    """Twist angle of bonds a1-a2 and a4-a3 around around the central bond a2-a3

    Args:
        a1,a2,a3,a4 (mdt.Atom): the atoms describing the dihedral
    """
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




