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
import math

import moldesign as mdt
from moldesign import units as u
from moldesign.external import transformations as trns

# NEWFEATURE: preliminary profiling indicates that, for interactive work, the UNITS LIBRARY is actually
#      the biggest bottleneck
# NEWFEATURE: pyramidalization aka out-of-plane bending

DIST_TOLERANCE = 1.0e-3 * u.angstrom
DIST_FORCE_CONSTANT = 1000.0 * u.kcalpermol / (u.angstrom**2)
ANGLE_TOLERANCE = 0.75 * u.degrees
ANGLE_FORCE_CONSTANT = 1500.0 * u.kcalpermol / (u.radians**2)

__all__ = 'distance set_distance angle set_angle dihedral set_dihedral'.split()

def distance(a1, a2):
    """ Return distance between two atoms

    Args:
        a1,a2 (mdt.Atom): the two atoms

    Returns:
        u.Scalar[length]: the distance
    """
    return a1.distance(a2)


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
    d = normalized(a1.position - a2.position)
    return d, -d


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
    costheta = np.dot(e12, e23)
    theta = safe_arccos(costheta)
    return theta * u.radians


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


def dihedral(a1, a2, a3, a4):
    """Twist angle of bonds a1-a2 and a4-a3 around around the central bond a2-a3

    Args:
        a1,a2,a3,a4 (mdt.Atom): the atoms describing the dihedral
    """
    r21 = (a1.position - a2.position).defunits_value()  # remove units immediately to improve speed
    r34 = (a4.position - a3.position).defunits_value()
    center_bond = (a2.position - a3.position).defunits_value()
    plane_normal = normalized(center_bond)
    r21_proj = r21 - plane_normal * r21.dot(plane_normal)
    r34_proj = r34 - plane_normal * r34.dot(plane_normal)
    va = normalized(r21_proj)
    vb = normalized(r34_proj)
    costheta = np.dot(va, vb)
    if np.allclose(costheta, 1.0):
        return 0.0 * u.radians
    elif np.allclose(costheta, -1.0):
        return u.pi * u.radians
    else:
        abstheta = safe_arccos(costheta)
        cross = np.cross(va, vb)
        theta = abstheta * np.sign(np.dot(cross, plane_normal))
        return (theta * u.radians) % (2.0 * u.pi * u.radians)


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
        assert a1.parent is not None
        assert a1.parent == a2.parent
    vec = a1.position - a2.position
    dist = np.sqrt(vec.dot(vec))
    direction = vec / dist
    delta = newlength - dist
    if np.abs(delta) < 1.0e-5 * delta.get_units(): return
    if not adjustmol:
        a1.position += direction * delta / 2.0
        a2.position -= direction * delta / 2.0
    else:
        mol = a1.parent
        indices, sign = _get_fragment_indices(mol, a1, a2)
        mol.positions[indices] = mol.positions[indices] + delta * direction * sign


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
        rotmat_l = trns.rotation_matrix(rotation / 2.0, axis, a2.position)
        rotmat_r = trns.rotation_matrix(-rotation / 2.0, axis, a2.position)

        a1.position = _apply_4trans(rotmat_l, a1.position)
        a3.position = _apply_4trans(rotmat_r, a3.position)

    else:
        mol = a2.parent
        indices, sign = _get_fragment_indices(mol, a1, a2)
        rotmat = trns.rotation_matrix(rotation, axis*sign, a2.position)
        mol.positions[indices] = _apply_4trans(rotmat, mol.positions[indices])


def set_dihedral(a1, a2, a3, a4, theta, adjustmol=True):
    """ Set the twist angle of atoms a1 and a4 around the central bond a2-a3. The atoms will be
    adjusted along the
    gradient of the angle.
    If ``adjustmol`` is True and the topology is unambiguous, then the entire molecule's positions
    will be modified as well

    Args:
        a1,a2,a3,a4 (mdt.Atom): atoms to adjust
        theta (u.Scalar[angle]): new angle to set
        adjustmol (bool): Adjust all atoms on either side of this bond?
    """
    # TODO: deal with co-linear a1/a4, a2, a3 - the angle is ill-defined \
    #      (should just an arbitrary axis normal to the central bond)
    current = dihedral(a1, a2, a3, a4)
    rotation = sub_angles(theta, current)
    if abs(rotation) < 1.0e-6: return

    axis = a2.position - a3.position
    if not adjustmol:
        rotmat_l = trns.rotation_matrix(-rotation / 2.0, axis, a3.position)
        rotmat_r = trns.rotation_matrix(rotation / 2.0, axis, a3.position)

        a1.position = _apply_4trans(rotmat_l, a1.position)
        a4.position = _apply_4trans(rotmat_r, a4.position)

    else:
        mol = a2.parent
        indices, sign = _get_fragment_indices(mol, a3, a2)
        rotmat = trns.rotation_matrix(rotation, axis*sign, a3.position)
        mol.positions[indices] = _apply_4trans(rotmat, mol.positions[indices])


class GeometryConstraint(object):
    """
    Base class - Keeps track of a 3D geometry constraint.
    The constraint is satisfied when self.current() == self.value
    """
    desc = 'base constraint'  # use this to identify the constraint in interfaces
    dof = 1  # number of degrees of freedom constrained (so that we can properly calculate temperature)

    def __init__(self, atoms, value=None, tolerance=None, force_constant=None):
        """ Initialization:

        Args:
            atoms (List[mdt.Atom]): atoms involved
            value (u.Scalar): constrain the coordinate to this value
            tolerance (u.Scalar): absolute tolerance (for iterative constraint enforcement)
            force_constant (u.Scalar[force]): optional, only for minimizations and/or use in
               restraints)
        """
        self.atoms = mdt.AtomList(atoms)
        self.mol = self.atoms[0].parent
        self.tolerance = tolerance
        self.force_constant = force_constant
        for atom in self.atoms:
            assert atom.parent is self.mol
        self.value = mdt.utils.if_not_none(value, self.current())

    def current(self):
        """
        Return the current value of the constrained quantity
        """
        raise NotImplementedError()

    def gradient(self):
        r"""
        Return the gradient of the constrained quantity
        Requires that self.atomgrad be implemented (or otherwise provided)

        .. math::
            \nabla G(\mathbf r)
        """
        grad = np.zeros(self.mol.ndims)
        atomvecs = self.atomgrad(*self.atoms)
        assert len(atomvecs) == len(self.atoms)
        grad = grad * atomvecs[0].get_units()
        for v, a in zip(atomvecs, self.atoms):
            grad[a.parent_slice] = v
        return grad

    def satisfied(self):
        """
        Returns:
            bool: True if this constraint is satisfied to within tolerance
        """
        return abs(self.error()) <= self.tolerance

    def restraint_penalty(self):
        """
        Returns:
            u.Scalar[energy]: energy penalty for restraints
        """
        return 0.5 * self.force_constant * self.error()**2

    def restraint_penalty_force(self):
        """
        Returns:
            u.Vector[energy]: forces from restraint
        """
        return -self.force_constant * self.gradient() * self.error()

    def error(self):
        """
        Returns:
            u.Scalar: deviation of current geometry from the constraint
        """
        return self.current() - self.value

    def __repr__(self):
        return '<{cls.__name__}{self.atoms}:value={self.value})>'.format(
            cls=self.__class__, self=self)

    def __str__(self):
        return 'Constraint: {self.desc}({atoms}) -> {self.value})>'.format(
            atoms=','.join([a.name for a in self.atoms]), self=self)


class DistanceConstraint(GeometryConstraint):
    desc = 'distance'
    atomgrad = staticmethod(distance_gradient)
    dof = 1

    def __init__(self, atom1, atom2, value=None,
                 tolerance=DIST_TOLERANCE, force_constant=DIST_FORCE_CONSTANT):
        self.a1 = atom1
        self.a2 = atom2
        super(DistanceConstraint, self).__init__([atom1, atom2], value=value,
                                                 tolerance=tolerance, force_constant=force_constant)

    def current(self):
        return self.a1.distance(self.a2)
    current.__doc__ = GeometryConstraint.current.__doc__


class AngleConstraint(GeometryConstraint):
    desc = 'angle'
    atomgrad = staticmethod(angle_gradient)
    dof = 1

    def __init__(self, atom1, atom2, atom3, value=None,
                 tolerance=ANGLE_TOLERANCE, force_constant=ANGLE_FORCE_CONSTANT):
        self.a1 = atom1
        self.a2 = atom2
        self.a3 = atom3
        super(AngleConstraint, self).__init__([atom1, atom2, atom3], value=value,
                                              tolerance=tolerance, force_constant=force_constant)

    def current(self):
        return angle(*self.atoms)
    current.__doc__ = GeometryConstraint.current.__doc__

    def error(self):
        return sub_angles(self.current(), self.value)
    error.__doc__ = GeometryConstraint.error.__doc__


class DihedralConstraint(GeometryConstraint):
    desc = 'dihedral'
    atomgrad = staticmethod(dihedral_gradient)
    dof = 1

    def __init__(self, atom1, atom2, atom3, atom4, value=None,
                 tolerance=ANGLE_TOLERANCE, force_constant=ANGLE_FORCE_CONSTANT):
        self.a1 = atom1
        self.a2 = atom2
        self.a3 = atom3
        self.a4 = atom4
        super(DihedralConstraint, self).__init__([atom1, atom2, atom3, atom4], value=value,
                                                 tolerance=tolerance, force_constant=force_constant)

    def current(self):
        return dihedral(*self.atoms)
    current.__doc__ = GeometryConstraint.current.__doc__

    def error(self):
        return sub_angles(self.current(), self.value)
    error.__doc__ = GeometryConstraint.error.__doc__


class FixedPosition(GeometryConstraint):
    """This constraint is different than the others because it's absolute, not relative"""
    desc = 'position'
    dof = 3

    def __init__(self, atom, point=None, value=None,
                 tolerance=DIST_TOLERANCE, force_constant=DIST_FORCE_CONSTANT):
        self.atom = atom
        self.point = mdt.utils.if_not_none(point, atom.position.copy())
        super(FixedPosition, self).__init__([atom], value=value,
                                            tolerance=tolerance, force_constant=force_constant)

    def current(self):
        return self.atom.position
    current.__doc__ = GeometryConstraint.current.__doc__

    def error(self):
        diff = self.atom.position - self.value
        return np.sqrt(diff.dot(diff))
    error.__doc__ = GeometryConstraint.error.__doc__

    def atomgrad(self):
        """
        Returns:
            u.Vector[length]: unit vector from the constraint location to the atom's actual
                location (len=3)
        """
        diff = self.atom.position - self.value
        return normalized(diff)


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
    result = mdt.AtomList(visited)
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

    indices = []
    for atom in frag:
        indices.append(range(atom.parent_slice.start, atom.parent_slice.stop))
    result = np.array(indices)
    _get_fragment_indices.cache[key] = (result, sign)
    return result, sign
_get_fragment_indices.cache = {}


def _apply_4trans(trans, vecs):
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
    mag = vec.dot(vec)
    if mag == 0.0: return vec * 0.0
    else: return vec / np.sqrt(vec.dot(vec))


def safe_arccos(costheta):
    """ Version of arccos that can handle numerical noise greater than 1.0
    """
    if abs(costheta) > 1.0:
        assert abs(costheta) - 1.0 < 1.0e-14
        return u.pi
    else:
        return np.arccos(costheta)


def sub_angles(a, b):
    """ Subtract two angles, keeping the result within [0,360)
    """
    c = a - b
    return (c + 180.0 * u.degrees) % (360.0 * u.degrees) - (180.0 * u.degrees)
