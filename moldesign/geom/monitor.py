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
        return self.CONSTRAINT(*self.atoms, **kwargs)


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
    COSTRAINT = constraints.AngleConstraint


@toplevel
class DihedralMonitor(Monitor):
    NUM_ATOMS = 4
    GETTER = staticmethod(coords.dihedral)
    SETTER = staticmethod(setcoord.set_dihedral)
    GRAD = staticmethod(grads.dihedral_gradient)
    CONSTRAINT = constraints.DihedralConstraint
