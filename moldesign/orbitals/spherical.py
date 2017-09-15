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

from ..mathutils import spherical_harmonics
from .. import units as u
from .. import utils
from . import Primitive, Gaussian, CartesianGaussian, PrimitiveSum


class SphericalGaussian(Primitive):
    r""" Stores a 3-dimensional real spherical gaussian function:

    .. math::
        G_{nlm}(\mathbf r) = C Y^l_m(\mathbf r - \mathbf r_0) r^l e^{-\alpha \left| \mathbf r - \mathbf r_0 \right|^2}

    where:
     - *C* is ``self.coeff``,
     - :math:`\alpha` is ``self.alpha``,
     - :math:`\mathbf r_0` is  ``self.center``,
     - _(l,m)_ are given by ``(self.l, self.m)``, and
     - :math:`Y^l_m(\mathbf{r})` are the orthonormal real-valued spherical harmonics.

    This function is evaluated as the the product of three terms that can be accessed separately:
      - ``self.coeff`` -  the constant prefactor _C_
      - ``self.radial_part(r)`` - :math:`e^{-a\left| \mathbf r - \mathbf r_0 \right|^2}`
      - ``self.angular_part(r)`` - :math:`Y^l_m(\mathbf r - \mathbf r_0) r^l`

    References:
        Schlegel and Frisch. Transformation between cartesian and pure spherical harmonic gaussians.
            Int J Quantum Chem 54, 83-87 (1995). doi:10.1002/qua.560540202

        https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
    """
    center = utils.Alias('radial_part.center')
    alpha = utils.Alias('radial_part.alpha')
    ndims = ndim = num_dims = 3

    def __init__(self, center, alpha, l, m, coeff=None, normalized=True):
        assert len(center) == 3, "Spherical gaussians must be 3-dimensional"
        self.l = l
        self.m = m
        self.radial_part = Gaussian(center, alpha, coeff=1.0, normalized=False)
        self._spherical_harmonic = spherical_harmonics.Y(l, m)
        super().__init__(coeff=coeff, normalized=normalized)

    def __repr__(self):
        return ("<3D Gaussian (Spherical) (coeff: {coeff:4.2f}, "
                "alpha: {alpha:4.2f}, "
                "(l,m) = {qnums}>").format(
                center=self.center, alpha=self.alpha, coeff=self.coeff,
                qnums=(self.l, self.m))

    def __mul__(self, other):
        raise NotImplementedError("Cannot yet multiply spherical gaussian functions")

    def _get_wfn_units(self):
        return u.MdtQuantity(self.coeff * u.default.length**self.l).units

    def __call__(self, coords):
        """ Evaluate this function at the given coordinates.

        Can be called either with a 1D column (e.g., ``[1,2,3]*u.angstrom ``) or
        an ARRAY of coordinates (``[[0,0,0],[1,1,1]] * u.angstrom``)

        Args:
            coords (Vector[length]): 3D Coordinates or list of 3D coordinates

        Returns:
            Scalar or Vector: function value(s) at the passed coordinates
        """
        result, disp, r2 = self.radial_part(coords, _getvals=True)
        result *= r2**(self.l/2.0) * self._spherical_harmonic(disp)
        return result * self.coeff

    def overlap(self, other, normalized=False):
        r""" Overlap of this function with another basis function

        .. math::
            \int f_1(\mathbf r) f_2(\mathbf r) d^N \mathbf r

        Args:
            other (Primitive):
            normalized (bool): If True, return the overlap of the two NORMALIZED functions.

        Returns:
            Scalar: value of the overlap
        """
        from scipy.special import factorial
        from ..data import ODD_FACTORIAL

        if (self.center != other.center).any():
            return self.to_cart().overlap(other)

        elif self.l != other.l or self.m != other.m:
            return 0.0

        else:
            newexp = self.alpha + other.alpha
            power = self.l + other.l + 2

            # In this case, radial part integrates to 1, so we just
            # integrate r^(2+l) exp(-(a1+a2) r^2) from 0 to infinity
            if power % 2 == 0:
                n = power // 2
                val = np.sqrt(np.pi/newexp) * ODD_FACTORIAL[2*n-1] / (newexp**n * 2**(n+1))
            else:
                n = (power - 1) // 2
                val = factorial(n, True) / (2 * newexp**(n+1))

        val *= self.coeff * other.coeff

        if normalized:
            return val / (self.norm * other.norm)
        else:
            return val

    def angular_part(self, coords):
        if len(coords.shape) > 1:
            axis = 1
        else:
            axis = None

        disp = coords - self.center
        prd = disp*disp
        r2 = prd.sum(axis=axis)
        return r2**(self.l/2.0) * self._spherical_harmonic(disp)

    def to_cart(self):
        """ Convert this function to a linear combination of cartesian functions

        Returns:
            List[CartesianGaussian]: cartesian components of this function
        """
        carts = [CartesianGaussian(self.center, self.alpha,
                                   cart_fn.powers, self.coeff * cart_fn.coeff,
                                   normalized=False)
                 for cart_fn in spherical_harmonics.SPHERE_TO_CART[self.l, self.m]]
        return PrimitiveSum(carts)

    @property
    def integral(self):
        return self.to_cart().integral
