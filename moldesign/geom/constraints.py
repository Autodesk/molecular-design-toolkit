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

import moldesign as mdt
import moldesign.molecules.atomcollections
from moldesign import units as u
from moldesign.mathutils import *
from .coords import *
from .grads import *

DIST_TOLERANCE = 1.0e-3 * u.angstrom
DIST_FORCE_CONSTANT = 1000.0 * u.kcalpermol / (u.angstrom**2)
ANGLE_TOLERANCE = 0.75 * u.degrees
ANGLE_FORCE_CONSTANT = 1500.0 * u.kcalpermol / (u.radians**2)


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
        self.atoms = moldesign.molecules.atomcollections.AtomList(atoms)
        self.mol = self.atoms[0].molecule
        self.tolerance = tolerance
        self.force_constant = force_constant
        for atom in self.atoms:
            assert atom.molecule is self.mol
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
            grad[a.index] = v
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


@toplevel
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


@toplevel
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


@toplevel
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


@toplevel
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