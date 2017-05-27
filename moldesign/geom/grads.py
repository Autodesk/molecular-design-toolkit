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
from moldesign.mathutils import normalized

from . import toplevel, angle


@toplevel
def distance_gradient(a1, a2):
    r""" Gradient of the distance between two atoms,

    .. math::
        \frac{\partial \mathbf{R}_1}{\partial \mathbf{r}} ||\mathbf{R}_1 - \mathbf{R}_2|| =
           \frac{\mathbf{R}_1 - \mathbf{R}_2}{||\mathbf{R}_1 - \mathbf{R}_2||}

    Args:
        a1,a2 (mdt.Atom): the two atoms

    Returns:
        Tuple[u.Vector[length], u.Vector[length]]: (gradient w.r.t. first atom, gradient w.r.t.
        second atom)
    """
    d = normalized(a1.position-a2.position)
    return d, -d


@toplevel
def angle_gradient(a1, a2, a3):
    r"""Gradient of the angle between bonds a2-a1 and a2-a3

    .. math::
        \nabla \theta_{ijkl} = \frac{\partial \theta_{ijkl}}{\partial \mathbf R}

    Args:
        a1,a2,a3 (mdt.Atom): the atoms describing the vector

    References:
        https://salilab.org/modeller/9v6/manual/node436.html
    """
    theta = angle(a1, a2, a3)
    costheta = np.cos(theta)
    p = np.power(1.0 - costheta**2, -0.5)
    vij = a1.position - a2.position
    vkj = a3.position - a2.position
    rij = np.sqrt(vij.dot(vij))
    rkj = np.sqrt(vkj.dot(vkj))
    eij = vij/rij
    ekj = vkj/rkj
    vec1 = p * (eij * costheta - ekj) / rij
    vec3 = p * (ekj * costheta - eij) / rkj
    vec2 = -vec1 - vec3
    return vec1, vec2, vec3


@toplevel
def dihedral_gradient(a1, a2, a3, a4):
    r""" Cartesian gradient of a dihedral coordinate,

    .. math::
        \nabla \theta_{ijkl} = \frac{\partial \theta_{ijkl}}{\partial \mathbf R}

    Args:
        a1,a2,a3,a4 (mdt.Atom): the atoms describing the dihedral

    References:
        https://salilab.org/modeller/9v6/manual/node436.html
    """
    vij = a1.position - a2.position
    vkj = a3.position - a2.position
    vkl = a3.position - a4.position
    vmj = vij.cross(vkj)
    vnk = vkj.cross(vkl)
    rkj = np.sqrt(vkj.dot(vkj))
    rmj = np.sqrt(vmj.dot(vmj))
    rnk = np.sqrt(vnk.dot(vnk))
    pijkj = vij.dot(vkj) / (rkj**2)
    pklkj = vkl.dot(vkj) / (rkj**2)

    vec1 = rkj * vmj / (rmj**2)
    vec4 = -rkj * vnk / (rnk**2)
    vec2 = vec1 * (pijkj - 1.0) - vec4 * pklkj
    vec3 = vec4 * (pklkj - 1.0) - vec1 * pijkj

    return -vec1 * u.radians, -vec2 * u.radians, -vec3 * u.radians, -vec4 * u.radians


def _atom_grad_to_mol_grad(atoms, grads):
    """ Convert list of gradients on atoms to a full-dimensional Nx3 gradient list (with 0s for
    uninvolved atoms)
    """
    m = atoms[0].molecule
    if len(grads) != len(atoms):
        raise ValueError('Number of gradients does not match number of atoms')
    mol_grad = np.zeros((m.num_atoms, 3))*grads[0].get_units()
    for v, a in zip(grads, atoms):
        mol_grad[a.index] = v
    return mol_grad
