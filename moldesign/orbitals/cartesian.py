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
import itertools

import numpy as np

from .. import units as u
from .. import utils
from . import Primitive, PrimitiveSum, Gaussian


class CartesianGaussian(Primitive):
    r""" Stores an N-dimensional gaussian function of the form:

    .. math::
        G(\mathbf r) = C \times \left( \prod_{i=1}^N{{r_i}^{p_i} } \right)
             e^{-\alpha |\mathbf r - \mathbf{r}_0|^2}

    For a three-dimensional gaussian, this is

    ..math::
        G(x,y,z) = C \times x^{p_1} y^{p_2} z^{p_3} e^{-\alpha |\mathbf r - \mathbf{r}_0|^2}

    where *C* is ``self.coeff``, :math:`\alpha` is ``self.alpha``, *r0* is ``self.center``, and
    :math:`p_1, p_2, ...` are given in the array ``self.powers``

    This function is evaluated as the the product of three terms that can be accessed separately:
      - ``self.coeff`` -  the constant prefactor _C_
      - ``self.radial_part(r)`` - :math:`e^{-\alpha \left| \mathbf r - \mathbf r_0 \right|^2}`
      - ``self.angular_part(r)`` - :math:`\prod_{i=1}^N{{r_i}^{p_i} }`

    Args:
        center (Vector[length]): location of the gaussian's centroid
        powers (List[int]): cartesian powers in each dimension (see
            equations in :class:`CartesianGaussian` docs)
        alpha (Scalar[1/length**2]): gaussian width parameter
        coeff (Scalar): multiplicative coefficient (if None, gaussian will be automatically
             normalized)

    Note:
        The dimensionality of the gaussian is determined by the dimensionality
        of the centroid location vector and the power vector. So, if scalars are passed for the
        ``center`` and ``powers``, it's 1-D. If length-3 vectors are passed for ``center``
        and ``powers``, it's 3D.

    References:
        Levine, Ira N. Quantum Chemistry, 5th ed. Prentice Hall, 2000. 486-94.
    """
    center = utils.Alias('radial_part.center')
    alpha = utils.Alias('radial_part.alpha')
    ndims = ndim = num_dims = utils.Alias('radial_part.ndim')

    def __init__(self, center, alpha, powers, coeff=None, normalized=True):
        assert len(powers) == len(center), "Inconsistent dimensionality - number of cartesian " \
                                           "powers must match dimensionality of centroid vector"

        self.powers = np.array(powers)
        self.shell = sum(self.powers)
        self.radial_part = Gaussian(center, alpha, coeff=1.0, normalized=False)
        super().__init__(coeff=coeff, normalized=normalized)

    def __repr__(self):
        return ("<{ndim}-D cartesian gaussian (norm: {norm:4.2f}, "
                "cartesian powers: {powers}, "
                "alpha: {exp:4.2f}, "
                "center: {center}>").format(
                ndim=self.ndim,
                center=self.center, exp=self.alpha,
                powers=tuple(self.powers), norm=self.norm)

    def _get_wfn_units(self):
        return u.MdtQuantity(self.coeff * u.default.length**self.angular).units

    @property
    def angular(self):
        """ Angular momentum of this function (sum of cartesian powers)
        """
        return self.powers.sum()

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
        olap = 0.0
        for other_prim in other.iterprimitives():
            if hasattr(other_prim, 'to_cart'):
                other_prim = other_prim.to_cart()

            for p in other_prim.iterprimitives():
                olap += self._overlap_cart_cart(p)

        if normalized:
            olap = olap / (self.norm * other.norm)

        return olap

    def _overlap_cart_cart(self, other):
        """  Overlap of two cartesian functions. Uses the PyQuante2 python implementation

        Users shouldn't need to call directly, this should be automatically called when necessary
        """
        from ..external.pyquante2.one import overlap
        lengthunits = u.default.get_default(self.center)
        alphaunits = (1/lengthunits**2).units

        alphas = [alphaunits.value_of(g.alpha) for g in (self, other)]
        centers = [lengthunits.value_of(g.center) for g in (self, other)]

        olap_base = overlap(alphas[0], self.powers, centers[0],
                       alphas[1], other.powers, centers[1])

        olap = self.coeff * other.coeff * olap_base

        if lengthunits != u.dimensionless:
            olap *= lengthunits ** (3 + sum(self.powers) + sum(other.powers))

        return olap

    def __call__(self, coords):
        """ Evaluate this function at the given coordinates.

        Can be called either with a 1D column (e.g., ``[1,2,3]*u.angstrom ``) or
        an ARRAY of coordinates (``[[0,0,0],[1,1,1]] * u.angstrom``)

        Args:
            coords (Vector[length, len=3] or Matrix[length, shape=(*,3)]): Coordinates or
                   list of 3D coordinates

        Examples:
            >>> g = CartesianGaussian([0,0,0]*u.angstrom, alpha=1.0/u.angstrom**2, powers=(0,0,0))
            >>> # Value at a single coordinate:
            >>> g([0,0,0] * u.angstrom)
            1.0
            >>> # Values at a list of coordinates
            >>> g[[0,0,0], [0,0,1], [0.5,0.5,0.5] * u.angstrom]
            array([ 1.0,  0.36787944,  0.47236655])

        Returns:
            Scalar or Vector: function value(s) at the passed coordinates
        """
        result = self.radial_part(coords)

        if self.shell > 0:
            result *= self.angular_part(coords)
        return result * self.coeff

    def angular_part(self, coords):
        if self.shell == 0:
            return 1.0

        if len(coords.shape) > 1:
            axis = 1
        else:
            axis = None

        disp = coords - self.center
        if hasattr(disp, 'units'):
            mag = disp.magnitude
            units = disp.units ** self.powers.sum()
        else:
            mag = disp
            units = 1.0
        return np.product(mag**self.powers, axis=axis) * units

    def __mul__(self, other):
        """ Returns product of two cartiesian gaussian functions as a list of product gaussians

        The returned gaussians all have the same center and width, and differ only in their
        angular parts.

        Args:
            other (CartesianGaussian): other gaussian basis function

        Returns:
            CartesianGaussian or PrimitiveSum[CartesianGaussian]: product functions
        """
        from scipy.special import binom

        if (self.center == other.center).all():  # this case is much easier than the general one
            cpy = self.copy()
            cpy.radial_part.alpha = self.alpha + other.alpha
            cpy.powers = self.powers + other.powers
            cpy.coeff = self.coeff * other.coeff
            return cpy

        newcenter = self.radial_part * other.radial_part

        partial_coeffs = [{} for idim in range(self.ndim)]
        for idim in range(self.ndim):
            r_self = self.center[idim]
            r_other = other.center[idim]
            r_new = newcenter.center[idim]

            for m in range(self.powers[idim]+1):
                for k in range(other.powers[idim]+1):
                    powercoeff = (binom(self.powers[idim], m) * binom(other.powers[idim], k) *
                                  ((r_new - r_self) ** (self.powers[idim]-m)) *
                                  ((r_new - r_other) ** (other.powers[idim]-k)))
                    if powercoeff == 0.0:
                        continue
                    newpower = m+k
                    if newpower not in partial_coeffs[idim]:
                        partial_coeffs[idim][newpower] = powercoeff
                    else:
                        partial_coeffs[idim][newpower] += powercoeff

        new_gaussians = []
        for powers in itertools.product(*[x.keys() for x in partial_coeffs]):
            final_coeff = self.coeff * other.coeff * newcenter.coeff
            for idim, p in enumerate(powers):
                final_coeff *= partial_coeffs[idim][p]
            new_gaussians.append(CartesianGaussian(newcenter.center, newcenter.alpha,
                                                   powers=powers, coeff=final_coeff,
                                                   normalized=False))

        if len(new_gaussians) == 0:  # no non-zero components, just return a zeroed gaussian
            return CartesianGaussian(newcenter, newcenter.alpha,
                                     self.powers + other.powers, coeff=0.0,
                                     normalized=False)
        elif len(new_gaussians) == 1:
            return new_gaussians[0]
        else:
            return PrimitiveSum(new_gaussians)

    @property
    def integral(self):
        r"""Integral of this this gaussian over all N-dimensional space.

        This is implemented only for 0 and positive integer cartesian powers.
        The integral is 0 if any of the powers are odd. Otherwise, the integral
        is given by:
        .. math::
            \int G(\mathbf r) d^N \mathbf r & = c \int d^N e^{-a x^2} \mathbf r
                    \prod_{i=1}^N{{r_i}^{p_i} }   \\
               &= (2a)^{-\sum_i p_i} \left( \frac{\pi}{2 a} \right) ^ {N/2} \prod_{i=1}^N{(p_i-1)!!}

        where *N* is the dimensionality of the gaussian, :math:`p_i` are the cartesian powers,
        and _!!_ is the "odd factorial" (:math:`n!!=1\times 3\times 5 \times ... \times n`)

        References:
            Dwight, Herbert B. Tables of Integrals and other Mathematical Data, 3rd ed.
                Macmillan 1957. 201.
        """
        from ..data import ODD_FACTORIAL

        integ = self.coeff * self.radial_part.integral
        for p in self.powers:
            if p == 0:  # no contribution
                continue
            elif p % 2 == 1:  # integral of odd function is exactly 0
                return 0.0
            elif p < 0:
                raise ValueError('Powers must be positive or 0')
            else:
                integ *= ODD_FACTORIAL[p-1]/((2.0 * self.alpha) ** (p//2))
        return integ
