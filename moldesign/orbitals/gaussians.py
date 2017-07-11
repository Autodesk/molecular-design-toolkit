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
import numpy as np

from .. import units as u


class AbstractFunction(object):
    """ Abstract base class for basis functions
    """

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
            other (AbstractFunction):
            normalized (bool): If True, return the overlap of the two NORMALIZED functions.

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

    def __str__(self):
        return "%d-D %s with norm %s" % (self.ndims, self.__class__, self.norm)


class Gaussian(AbstractFunction):
    r""" N-dimensional gaussian function.

    The function is given by:
    .. math::
        G(\mathbf r) = C e^{-a\left| \mathbf r - \mathbf r_0 \right|^2}

    where *C* is ``self.coeff``, *a* is ``self.exp``, and :math:`\mathbf r_0` is ``self.center``.
    """
    def __init__(self, center, exp, coeff=None):
        self.center = u.array(center)
        self.exp = exp

        if coeff is None:
            self.coeff = 1.0  # dummy value overwritten by self.normalize()
            self.normalize()
        else:
            self.coeff = coeff

    def __call__(self, coords, _getvals=False):
        """ Evaluate this function at the given coordinates.

        Can be called either with a 1D column (e.g., ``[1,2,3]*u.angstrom ``) or
        an ARRAY of coordinates (``[[0,0,0],[1,1,1]] * u.angstrom``)

        Args:
            _include_angular (bool): include the contribution from the non-exponential parts
                (for computational efficiency, this can be omitted now and included later)

        Args:
            coords (Vector[length]): 3D Coordinates or list of 3D coordinates

        Returns:
            Scalar or Vector: function value(s) at the passed coordinates
        """
        if len(coords.shape) > 1:
            axis = 1
        else:
            axis = None

        disp = coords - self.center
        prd = disp*disp  # don't use np.dot - allow multiple coords at once
        r2 = prd.sum(axis=axis)

        result = self.coeff * np.exp(-self.exp * r2)
        if _getvals:
            return result, disp, r2  # extra quantities so that they don't need to be recomputed
        else:
            return result

    def __repr__(self):
        return ("<{ndim}-D gaussian (coeff: {coeff:4.2f}, "
                "exponent: {exp:4.2f}, "
                "center: {center}>").format(
                ndim=self.ndim,
                center=self.center, exp=self.exp,
                coeff=self.coeff)

    @property
    def ndim(self):
        return len(self.center)
    num_dimensions = ndims = ndim

    def __mul__(self, other):
        """ Returns a new gaussian centroid.
        """
        a = self.exp
        b = other.exp
        exp = a + b
        center = (a*self.center + b*other.center)/(a+b)
        prefactor = self.coeff * other.coeff
        for i in range(self.ndim):
            prefactor *= np.exp(-(a*self.center[i]**2 + b*other.center[i]**2) +
                                exp * center[i]**2)

        return Gaussian(center, exp, coeff=prefactor)

    @property
    def integral(self):
        r"""Integral of this this gaussian over all N-dimensional space
        """
        return self.coeff * (np.pi/self.exp)**(self.ndim/2.0)


def cart_to_powers(s):
    """ Convert a string to a list of cartesian powers

    Examples:
        >>> cart_to_powers('y')
        [0, 1, 0]
        >>> cart_to_powers('xxyz')
        [2, 1, 1]
        >>> cart_to_powers('zx^3')
        [3,0,1]
    """
    powers = [0, 0, 0]
    chars = iter(s)
    lastchar = None
    while True:
        try:
            char = next(chars)
        except StopIteration:
            break

        if char == '^':
            power = int(next(chars))
            powers[DIMLABELS[lastchar]] += power - 1
            lastchar = None
        else:
            powers[DIMLABELS[char]] += 1
            lastchar = char

    return powers


DIMLABELS = {'x': 0, 'y': 1, 'z': 2}
