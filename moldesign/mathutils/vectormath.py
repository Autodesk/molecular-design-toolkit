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

from .. import units as u
from ..utils import exports


@exports
def perpendicular(vec):
    """ Return arbitrary unit vector(s) perpendicular to one or more vectors.

    This obviously doesn't have a unique solution, but is useful for various computations.

    Args:
        vec (Vector or List[Vector]): 3-D vector or list thereof

    Returns:
        vec (np.ndarray): normalized unit vector in a perpendicular direction
    """
    vectorized = len(vec.shape) > 1
    direction = normalized(vec)
    if vectorized:
        cross_axis = np.array([[0.0, 0.0, 1.0] if d[2] < 0.9 else [0.0, 1.0, 0.0]
                              for d in direction])
    else:
        if abs(direction[2]) < 0.9:
            cross_axis = np.array([0.0, 0.0, 1.0])
        else:
            cross_axis = np.array([0.0, 1.0, 0.0])
    perp = normalized(np.cross(direction, cross_axis))
    return perp


@exports
def norm(vec):
    """ Calculate norm of a vector or list thereof

    Args:
        vec (Vector or List[Vector]): vector(s) to compute norm of

    Returns:
        Scalar or List[Scalar]: norms
    """
    if len(vec.shape) == 1:  # it's just a single column vector
        return np.sqrt(vec.dot(vec))
    else:  # treat as list of vectors
        return np.sqrt((vec*vec).sum(axis=1))


@exports
def normalized(vector, zero_as_zero=True):
    """ Create normalized versions of a vector or lists of vectors.

    Args:
        vector (Vector or List[Vector]): vector(s) to be normalized
        zero_as_zero (bool): if True, return a 0-vector if a 0-vector is passed; otherwise, will
           follow default system behavior (depending on numpy's configuration)

    Returns:
        Vector or List[Vector]: normalized vector(s)
    """
    vec = getattr(vector, 'magnitude', vector)  # strip units right away if necessary
    mag = norm(vec)
    if len(vec.shape) == 1:  # it's just a single column vector
        if mag == 0.0 and zero_as_zero:
            return vec*0.0
        else:
            return vec/mag
    else:  # treat as list of vectors
        if zero_as_zero:
            mag[mag == 0.0] = 1.0  # prevent div by 0
        return vec / mag[:, None]


@exports
def alignment_rotation(v1, v2, handle_linear=True):
    """ Calculate rotation angle(s) and axi(e)s to make v1 parallel with v2

    Args:
        v1 (vector or List[Vector]): 3-dimensional vector(s) to create rotation for
        v2 (vector or List[Vector]): 3-dimensional vector(s) to make v1 parallel to
        handle_linear (bool): if v1 is parallel or anti-parallel to v2, return an arbitrary
           vector perpendicular to both as the axis (otherwise, returns a 0-vector)

    Returns:
        MdtQuantity[angle]: angle between the two
        np.ndarray[len=3]: rotation axis (unit vector)

    References:
        https://stackoverflow.com/a/10145056/1958900
    """
    e1 = normalized(v1)
    e2 = normalized(v2)

    vectorize = len(e1.shape) > 1

    normal = np.cross(e1, e2)
    s = norm(normal)
    if vectorize:
        c = (e1*e2).sum(axis=1)
        if handle_linear:
            linear_indices = s == 0.0
            if linear_indices.any():
                normal[linear_indices] = perpendicular(v1[linear_indices])
                s[linear_indices] = 1.0
    else:
        c = np.dot(e1, e2)
        if handle_linear and s == 0.0:
            normal = perpendicular(v1)
            s = 1.0

    angle = np.arctan2(s, c)
    return angle*u.radian, normal / s


@exports
def safe_arccos(costheta):
    """ Version of arccos that can handle numerical noise greater than 1.0
    """
    if hasattr(costheta, 'shape') and costheta.shape:  # vector version
        assert (np.abs(costheta)-1.0 < 1.0e-13).all()
        costheta[costheta > 1.0] = 1.0
        costheta[costheta < -1.0] = -1.0
        return np.arccos(costheta)

    else:
        if abs(costheta) > 1.0:
            assert abs(costheta) - 1.0 < 1.0e-14
            return u.pi
        else:
            return np.arccos(costheta)


@exports
def sub_angles(a, b):
    """ Subtract two angles, keeping the result within [-180,180)
    """
    return normalize_angle(a - b)


@exports
def normalize_angle(c):
    """ Normalize an angle to the interval [-180,180)
    """
    return (c + 180.0 * u.degrees) % (360.0 * u.degrees) - (180.0 * u.degrees)


@exports
def apply_4x4_transform(trans, vecs):
    """
    Applies a 4x4 transformation vector so one or more 3-D position vector
    :param trans:
    :param vecs:
    :return: transformed position vector
    """
    has_units = False
    if hasattr(vecs, 'get_units'):
        has_units = True
        units = vecs.get_units()
        vecs = vecs.magnitude
    if len(vecs.shape) == 1:
        v = np.ones(4)
        v[:3] = vecs
        vt = trans.dot(v)
        result = vt[:3] / vt[3]
    else:
        v = np.ones((4, len(vecs)))
        v[:3, :] = vecs.T
        vt = trans.dot(v)
        result = (vt[:3] / vt[3]).T
    if has_units:
        result = result * units
    return result
