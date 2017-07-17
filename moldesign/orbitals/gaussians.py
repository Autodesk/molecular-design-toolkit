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

from .. import units as u
from . import Primitive


class Gaussian(Primitive):
    r""" N-dimensional gaussian function.

    The function is given by:
    .. math::
        G(\mathbf r) = C e^{-\alpha\left| \mathbf r - \mathbf r_0 \right|^2}

    where *C* is ``self.coeff``, *a* is ``self.alpha``, and :math:`\mathbf r_0` is ``self.center``.

    Args:
        center (Vector): centroid of this function (r_0 in the equation above)
        alpha (Scalar): constant in the exponential (alpha in the equation above)
        coeff (Scalar): coefficient for this function. If normalized=True, the coefficient
           multiplies the _normalized_ gaussian. If not passed, the gaussian will automatically
           be normalized
        normalized (bool): if True, normalize the function before applying the coefficient
    """
    def __init__(self, center, alpha, coeff=None, normalized=True):
        self.center = u.default.convert(center)
        self.alpha = u.default.convert(alpha)
        super().__init__(coeff=coeff, normalized=normalized)

    def _get_wfn_units(self):
        return u.MdtQuantity(self.coeff).units

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

        result = self.coeff * np.exp(-self.alpha * r2)
        if _getvals:
            return result, disp, r2  # extra quantities so that they don't need to be recomputed
        else:
            return result

    def __repr__(self):
        return ("<{ndim}-D gaussian (coeff: {coeff:4.2f}, "
                "exponent: {exp:4.2f}, "
                "center: {center}>").format(
                ndim=self.ndim,
                center=self.center, exp=self.alpha,
                coeff=self.coeff)

    @property
    def ndim(self):
        return len(self.center)
    num_dimensions = ndims = ndim

    def __mul__(self, other):
        """ Returns a new gaussian centroid.
        """
        a1 = self.alpha
        a2 = other.alpha
        x1 = self.center
        x2 = other.center
        alpha = a1 + a2
        center = (a1*x1 + a2*x2)/(a1+a2)
        prefactor = self.coeff * other.coeff
        for i in range(self.ndim):
            prefactor *= np.exp(-(a1*x1[i]**2 + a2*x2[i]**2) +
                                alpha * center[i]**2)

        return Gaussian(center, alpha, coeff=prefactor, normalized=False)

    @property
    def integral(self):
        r"""Integral of this this gaussian over all N-dimensional space
        """
        return self.coeff * (np.pi/self.alpha)**(self.ndim/2.0)

    def to_cart(self):
        from . import CartesianGaussian
        return CartesianGaussian(self.center, self.alpha, (0,0,0),
                                 self.coeff, normalized=False)


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
