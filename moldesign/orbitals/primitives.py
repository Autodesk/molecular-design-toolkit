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
import copy
import itertools

import numpy as np
from .. import utils


class Primitive(object):
    """ Abstract base class for basis function primitives. All functions are assumed real.
    """
    def __init__(self, coeff=None, normalized=False):
        self.coeff = 1.0  # initialize with a dummy value

        if normalized or (coeff is None):
            self.normalize()

        if coeff is not None:
            self.coeff *= coeff

    def _get_wfn_units(self):
        return NotImplementedError()

    def copy(self):
        return copy.deepcopy(self)

    def normalize(self):
        """ Give this function unit norm by adjusting its coefficient
        """
        self.coeff /= self.norm

    def overlap(self, other, normalized=False):
        r""" Overlap of this function with another:

        .. math::
            \int f_1(\mathbf r) f_2(\mathbf r) d^N \mathbf r

        Args:
            other (Primitive):
            normalized (bool): If True, return the overlap of the two NORMALIZED functions.

        Note:
            This implementation computes overlap by computing the integral of the product of the
            two functions. Subclasses may override this method if more efficient techniques
            are applicable

        Returns:
            Scalar: value of the overlap
        """
        newfn = self * other
        integral = newfn.integral
        if normalized:
            integral /= (self.norm*other.norm)
        return integral

    @property
    def norm(self):
        r""" The L2-Norm of this object, calculated as the square root of its self overlap.

        .. math::
            \sqrt{\int \left| G(\mathbf r) \right|^2 d^N \mathbf r}
        """
        return np.sqrt(self.overlap(self))

    def iterprimitives(self):
        yield self


class PrimitiveSum(object):
    """ Stores a linear combination of primitive functions.

    Args:
        primitives (List[Primitives]): List of primitives, if available
    """
    __getitem__ = utils.Alias('primitives.__getitem__')
    __len__ = utils.Alias('primitives.__len__')
    __iter__ = utils.Alias('primitives.__iter__')
    iterprimitives = __iter__

    def __init__(self, primitives):
        self.primitives = primitives
        self._wfn_units = self.primitives[0]._get_wfn_units()

    def _get_wfn_units(self):
        return self._wfn_units

    def __call__(self, coords):
        outvals = np.zeros(len(coords)) * self._get_wfn_units()
        for primitive in self.primitives:
            outvals += primitive(coords)
        return outvals

    @property
    def num_primitives(self):
        return len(self.primitives)

    @property
    def center(self):
        c = self.primitives[0].center.copy()
        for p in self:
            if (p.center != c).any():
                raise ValueError("This is a sum over multiple centers")
        else:
            return c

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
            primitive.coeff *= prefactor

    def overlap(self, other):
        """ Calculate orbital overlap with another object

        Args:
            other (AbstractFunction or Orbital): object to calculate overlaps with
        """
        olap = sum(p1.overlap(p2)
                   for p1, p2 in itertools.product(self.iterprimitives(),
                                                   other.iterprimitives()))
        return olap

    @property
    def integral(self):
        return sum(x.integral for x in self)

    def __str__(self):
        return '%s with %d primitives' % (self.__class__.__name__, self.num_primitives)

    def __repr__(self):
        return '<%s>' % self

    def __mul__(self, other):
        newprims = [p1*p2
                    for p1, p2 in itertools.product(self.iterprimitives(), other.iterprimitives())]
        return PrimitiveSum(newprims)
