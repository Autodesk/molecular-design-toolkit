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

from moldesign import units as u
from moldesign.utils import exports


def perpendicular(vec):
    """ Return an arbitrary unit vector perpendicular to vec.

    This obviously doesn't have a unique solution, but is useful for various computations.

    Args:
        vec (Vector): 3-D vector

    Returns:
        vec (np.ndarray): normalized unit vector in a perpendicular direction
    """

    assert vec.shape == (3,)
    direction = normalized(vec)
    if abs(direction[2]) < 0.9:
        cross_axis = np.array([0.0, 0.0, 1.0])
    else:
        cross_axis = np.array([0.0, 1.0, 0.0])
    perp = normalized(np.cross(direction, cross_axis))
    return perp


def norm(vec):
    """ Calculate norm of a vector

    Args:
        vec (Vector or List[Vector]): vector(s) to compute norm of

    Returns:
        Scalar or List[Scalar]: norms
    """
    if len(vec.shape) == 1:  # it's just a single column vector
        return np.sqrt(vec.dot(vec))
    else:  # treat as list of vectors
        return np.sqrt((vec*vec).sum(axis=1))


def normalized(vec, zero_as_zero=True):
    """ Create normalized versions of a vector or lists of vectors.

    Args:
        vec (Vector or List[Vector]): vector(s) to be normalized
        zero_as_zero (bool): if True, return a 0-vector if a 0-vector is passed; otherwise, will
           follow default system behavior (depending on numpy's configuration)

    Returns:
        Vector or List[Vector]: normalized vector(s)
    """
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


def alignment_rotation(v1, v2):
    """ Calculate rotation angle and axis to make v1 parallel with v2

    Args:
        v1 (Vector): 3-dimensional vector to create rotation for
        v2 (Vector): 3-dimensional vector to make v1 parallel to

    Returns:
        MdtQuantity[angle]: angle between the two
        np.ndarray[len=3]: rotation axis (unit vector)

    References:
        https://stackoverflow.com/a/10145056/1958900
    """
    e1 = normalized(v1)
    e2 = normalized(v2)

    normal = np.cross(e1, e2)
    s = norm(normal)
    c = np.dot(e1, e2)
    angle = np.arctan2(s, c)
    return angle*u.radian, normalized(normal)


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


def sub_angles(a, b):
    """ Subtract two angles, keeping the result within [-180,180)
    """
    c = a - b
    return (c + 180.0 * u.degrees) % (360.0 * u.degrees) - (180.0 * u.degrees)


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