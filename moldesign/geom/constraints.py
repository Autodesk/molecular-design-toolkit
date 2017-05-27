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

import moldesign as mdt
import moldesign.molecules.atomcollections
from moldesign import units as u
from moldesign.mathutils import *
from .coords import *
from .grads import *
from .grads import _atom_grad_to_mol_grad

DIST_TOLERANCE = 1.0e-6 * u.angstrom
DIST_FORCE_CONSTANT = 1000.0 * u.kcalpermol / (u.angstrom**2)
ANGLE_TOLERANCE = 1.0e-4 * u.degrees
ANGLE_FORCE_CONSTANT = 1500.0 * u.kcalpermol / (u.radians**2)


class GeometryConstraint(object):
    """ Base class - Keeps track of a 3D geometry constraint.
    The constraint is satisfied when self.current() == self.value

    Args:
        atoms (List[mdt.Atom]): atoms involved
        value (u.Scalar): constrain the coordinate to this value
        tolerance (u.Scalar): absolute tolerance (for iterative constraint enforcement)
        force_constant (u.Scalar[force]): optional, only for minimizations and/or use in
           restraints)
    """
    desc = 'base constraint'  # use this to identify the constraint in interfaces
    dof = 1  # number of degrees of freedom constrained (so that we can properly calculate temperature)

    def __init__(self, atoms, value=None, tolerance=None, force_constant=None):
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
        return _atom_grad_to_mol_grad(self.atoms, self.atomgrad(*self.atoms))

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

    def _constraintsig(self):
        """ Returns a unique key that lets us figure out if we have duplicate or conflicting
        constraints
        """
        return tuple([self.desc] + [atom.index for atom in self.atoms])


class DistanceConstraint(GeometryConstraint):
    desc = 'distance'
    atomgrad = staticmethod(distance_gradient)
    dof = 1

    def __init__(self, atom1, atom2, value=None,
                 tolerance=DIST_TOLERANCE, force_constant=DIST_FORCE_CONSTANT):
        self.a1 = atom1
        self.a2 = atom2
        super().__init__([atom1, atom2], value=value,
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
        super().__init__([atom1, atom2, atom3], value=value,
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
        super().__init__([atom1, atom2, atom3, atom4], value=value,
                                                 tolerance=tolerance, force_constant=force_constant)

    def current(self):
        return dihedral(*self.atoms)
    current.__doc__ = GeometryConstraint.current.__doc__

    def error(self):
        return sub_angles(self.current(), self.value)
    error.__doc__ = GeometryConstraint.error.__doc__


class FixedPosition(GeometryConstraint):
    """Fixes a single atom at a given location

    Note:
        The gradient of this function is singular and discontinuous when the constraint is satisfied,
        leading to poor results in iterative methods such as SHAKE.

        In such cases, this constraint should be automatically replaced with three
        :class:`FixedCoordinate` constraints on the atom's x, y, and z coordinates.
    """
    desc = 'position'
    dof = 3

    def __init__(self, atom, value=None,
                 tolerance=DIST_TOLERANCE, force_constant=DIST_FORCE_CONSTANT):
        self.atom = atom
        if value is None:
            self.value = mdt.utils.if_not_none(value, atom.position.copy())
        else:
            self.value = value.copy()
        super().__init__([atom], value=self.value,
                                            tolerance=tolerance, force_constant=force_constant)

    def current(self):
        return self.atom.position
    current.__doc__ = GeometryConstraint.current.__doc__

    def error(self):
        diff = self.atom.position - self.value
        return np.sqrt(diff.dot(diff))
    error.__doc__ = GeometryConstraint.error.__doc__

    def atomgrad(self, atom=None):
        """
        Note:
            For numerical reasons, this returns a vector in the [1,1,1] direction if the
            constraint is exactly met

        Returns:
            u.Vector[length]: unit vector from the constraint location to the atom's actual
                location (len=3)
        """
        if atom: assert atom is self.atom
        diff = self.atom.position - self.value
        grad = normalized(diff)
        if (grad.magnitude == np.zeros(3)).all():
            grad.magnitude[:] = np.ones(3) / np.sqrt(3)
        return [grad]


class FixedCoordinate(GeometryConstraint):
    """Fixes a single, linear degree of freedom for an atom, so that it's constrained to a plane.

    Args:
        atom (moldesign.Atom): atom to constrain
        vector (np.ndarray): direction to constrain
        value (units.Scalar[length]): constraint value (i.e., we constrain the dot product of the
            atom's position and the normalized direction vector to be this value)
    """
    desc = 'coordinate'
    dof = 1

    def __init__(self, atom, vector, value=None,
                 tolerance=DIST_TOLERANCE, force_constant=DIST_FORCE_CONSTANT):
        self.atom = atom
        self.vector = normalized(vector)
        if value is None:
            self.value = mdt.utils.if_not_none(value, self.current())
        else:
            self.value = value.copy()
        super().__init__([atom], value=self.value,
                                              tolerance=tolerance, force_constant=force_constant)

    def current(self):
        return self.atom.position.dot(self.vector)
    current.__doc__ = GeometryConstraint.current.__doc__

    def atomgrad(self, atom=None):
        return [self.vector] * u.ureg.dimensionless  # that was easy

    def _constraintsig(self):
        return super()._constraintsig() + (self.vector,)
