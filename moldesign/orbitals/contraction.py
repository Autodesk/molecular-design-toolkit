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
from .orbitals import Orbital


class BasisContraction(Orbital):
    """ Stores a linear combination of primitive functions.

    Args:
        primitives (List[PrimitiveBase]): List of primitives, if available
    """
    def __init__(self, primitives):
        self.primitives = primitives

    def __call__(self, coords):
        outvals = np.zeros(len(coords))
        for primitive in self.primitives:
            outvals += primitive(coords)
        return outvals

    @property
    def num_primitives(self):
        return len(self.primitives)

    @property
    def norm(self):
        """ Scalar: :math:`\sqrt{<i|i>}`
        """
        norm = 0.0
        for p1 in self.primitives:
            for p2 in self.primitives:
                norm += p1.overlap(p2)
        return np.sqrt(norm)

    def normalize(self):
        """ Scale primitive coefficients to normalize this basis function
        """
        prefactor = 1.0 / self.norm
        for primitive in self.primitives:
            primitive *= prefactor

    def overlap(self, other):
        """ Calculate orbital overlap with another object

        Args:
            other (AbstractFunction or Orbital): object to calculate overlaps with
        """
        olap = sum(p1.overlap(other) for p1 in self.primitives)
        return olap

    def __str__(self):
        return '%s with %d primitives' % (self.__class__.__name__, self.num_primitives)

    def __repr__(self):
        return '<%s>' % self