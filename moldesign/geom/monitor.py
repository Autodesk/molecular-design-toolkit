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

from . import toplevel
from . import constraints, grads, coords, setcoord


class Monitor(object):
    def __init__(self, *atoms):
        if len(atoms) != self.NUM_ATOMS:
            raise ValueError('%s requires %d atoms, but %d passed' %
                             (type(self), self.NUM_ATOMS, len(atoms)))
        self.atoms = atoms

    @property
    def value(self):
        return self.GETTER(*self.atoms)

    @value.setter
    def value(self, val):
        args = self.atoms + (val,)
        self.SETTER(*args)

    def gradient(self):
        return grads._atom_grad_to_mol_grad(self.atoms, self.GRAD(*self.atoms))

    @mdt.utils.kwargs_from(constraints.GeometryConstraint)
    def constrain(self, **kwargs):
        """ Constrain this coordinate.

        This will add a new item to the parent molecule's constraint list.

        Args:
            **kwargs (dict): kwargs for constraints.GeometryConstraint

        Returns:
            constraints.GeometryConstraint: the constraint object
        """
        c = self.CONSTRAINT(*self.atoms, **kwargs)
        mol = self.atoms[0].molecule
        for atom in mol.atoms[1:]:
            if atom.molecule is not mol:
                raise ValueError("Can't create constraint; atoms are not part of the same Molecule")
        mol.constraints.append(c)
        mol._reset_methods()
        return c

    def __call__(self, obj):
        """ Calculate this value for the given trajectory

        Args:
            obj (mdt.Molecule or mdt.Trajectory): molecule or trajectory to measure

        Returns:
            moldesign.units.Quantity: this coordinate's value (for a molecule), or a list of values
               (for a trajectory)

        Note:
            Atoms are identified by their index only; the atoms defined in the Monitor must have
            the same indices as those in the passed object
        """
        return self.GETTER(*(obj.atoms[a.index] for a in self.atoms))

    def __str__(self):
        return '%s: %s' % (type(self).__name__, self.value)

    def __repr__(self):
        return '<%s for atoms %s: %s>' % (type(self).__name__,
                                                ','.join(str(atom.index) for atom in self.atoms),
                                                self.value)


@toplevel
class DistanceMonitor(Monitor):
    NUM_ATOMS = 2
    GETTER = staticmethod(coords.distance)
    SETTER = staticmethod(setcoord.set_distance)
    GRAD = staticmethod(grads.distance_gradient)
    CONSTRAINT = constraints.DistanceConstraint


@toplevel
class AngleMonitor(Monitor):
    NUM_ATOMS = 3
    GETTER = staticmethod(coords.angle)
    SETTER = staticmethod(setcoord.set_angle)
    GRAD = staticmethod(grads.angle_gradient)
    CONSTRAINT = constraints.AngleConstraint


@toplevel
class DihedralMonitor(Monitor):
    def __init__(self, *atoms):
        if len(atoms) in (1, 2):
            atoms = coords._infer_dihedral(*atoms)
        super(DihedralMonitor, self).__init__(*atoms)

    NUM_ATOMS = 4
    GETTER = staticmethod(coords.dihedral)
    SETTER = staticmethod(setcoord.set_dihedral)
    GRAD = staticmethod(grads.dihedral_gradient)
    CONSTRAINT = constraints.DihedralConstraint
