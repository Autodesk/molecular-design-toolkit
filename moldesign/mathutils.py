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


def perpendicular(vec):
    assert vec.shape == (3,)
    direction = normalized(vec)
    if abs(direction[2]) < 0.9:
        cross_axis = np.array([0.0, 0.0, 1.0])
    else:
        cross_axis = np.array([0.0, 1.0, 0.0])
    perp = normalized(np.cross(direction, cross_axis))
    return perp


def normalized(vec):
    """ Return a vector normalized in L2.
    If vector is 0, return 0

    Args:
        vec (u.Vector): vector to be normalized

    Returns:
        u.Vector: normalized vector
    """
    if len(vec.shape) == 1:  # it's just a single column vector
        mag = vec.dot(vec)
        if mag == 0.0:
            return vec*0.0
        else:
            return vec/np.sqrt(mag)
    else:  # treat as list of vectors
        mag = (vec*vec).sum(axis=1)
        mag[mag == 0.0] = 1.0  # prevent div by 0
        return vec / np.sqrt(mag)[:, None]

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