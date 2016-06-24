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

    return vec1 * u.radians, vec2 * u.radians, vec3 * u.radians, vec4 * u.radians