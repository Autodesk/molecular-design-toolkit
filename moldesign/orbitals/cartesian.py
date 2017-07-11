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
from scipy.special import binom

from .. import units as u
from . import Gaussian


class CartesianGaussian(Gaussian):
    r""" Stores an N-dimensional gaussian function of the form:

    .. math::
        G(\mathbf r) = C \times \left( \prod_{i=1}^N{{r_i}^{p_i} } \right)
             e^{-a |\mathbf r - \mathbf{r}_0|^2}

    For a three-dimensional gaussian, this is

    ..math::
        G(x,y,z) = C \times x^{p_1} y^{p_2} z^{p_3} e^{-a |\mathbf r - \mathbf{r}_0|^2}

    where *C* is ``self.coeff``, *a* is ``self.exp``, *r0* is ``self.center``, and
    :math:`p_1, p_2, ...` are given in the array ``self.powers``

    References:
        Levine, Ira N. Quantum Chemistry, 5th ed. Prentice Hall, 2000. 486-94.

    Args:
        center (Vector[length]): location of the gaussian's centroid
        powers (List[int]): cartesian powers in each dimension (see
            equations in :class:`CartesianGaussian` docs)
        exp (Scalar[1/length**2]): gaussian width parameter
        coeff (Scalar): multiplicative coefficient (if None, gaussian will be automatically
             normalized)

    Note:
        The dimensionality of the gaussian is determined by the dimensionality
        of the centroid location vector and the power vector. So, if scalars are passed for the
        ``center`` and ``powers``, it's 1-D. If length-3 vectors are passed for ``center``
        and ``powers``, it's 3D.
    """
    def __init__(self, center, exp, powers, coeff=None):
        assert len(powers) == len(center), "Inconsistent dimensionality - number of cartesian " \
                                           "powers must match dimensionality of centroid vector"

        super().__init__(center, exp, coeff=1.0)  # temporarily set coefficient
        self.powers = np.array(powers)
        self.shell = sum(self.powers)
        if coeff is None:
            self.normalize()
        else:
            self.coeff = coeff

    def __repr__(self):
        return ("<{ndim}-D cartesian gaussian (coeff: {coeff:4.2f}, "
                "cartesian powers: {powers}, "
                "exponent: {exp:4.2f}, "
                "center: {center}>").format(
                ndim=self.ndim,
                center=self.center, exp=self.exp,
                powers=tuple(self.powers), coeff=self.coeff)

    @property
    def angular(self):
        """ Angular momentum of this function (sum of cartesian powers)
        """
        return self.powers.sum()

    def __call__(self, coords, _include_angular=True):
        """ Evaluate this function at the given coordinates.

        Can be called either with a 1D column (e.g., ``[1,2,3]*u.angstrom ``) or
        an ARRAY of coordinates (``[[0,0,0],[1,1,1]] * u.angstrom``)

        Args:
            coords (Vector[length, len=3] or Matrix[length, shape=(*,3)]): Coordinates or
                   list of 3D coordinates
            _include_cartesian (bool): include the contribution from the cartesian components
                (for computational efficiency, this can sometimes omited now and included later)

        Examples:
            >>> g = CartesianGaussian([0,0,0]*u.angstrom, exp=1.0/u.angstrom**2, powers=(0,0,0))
            >>> # Value at a single coordinate:
            >>> g([0,0,0] * u.angstrom)
            1.0
            >>> # Values at a list of coordinates
            >>> g[[0,0,0], [0,0,1], [0.5,0.5,0.5] * u.angstrom]
            array([ 1.0,  0.36787944,  0.47236655])

        Returns:
            Scalar or Vector: function value(s) at the passed coordinates
        """
        result = super().__call__(coords)

        if self.shell > 0 and _include_angular:
            result *= self.angular_part(coords)
        return result

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
            other (CartesianGaussian): other gaussian wavepacket

        Returns:
            CartesianGaussian or List[CartesianGaussian]: product gaussian(s)

        TODO:
           - Should return a composite gaussian object with the same interface as regular gaussians,
             rather than either an object OR a list
        """

        if (self.center == other.center).all():  # this case is much easier than the general one
            cpy = self.copy()
            cpy.exp = self.exp + other.exp
            cpy.powers = self.powers + other.powers
            cpy.coeff = self.coeff * other.coeff
            return cpy

        newcenter = super().__mul__(other)

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
            final_coeff = 1.0 * newcenter.coeff
            for idim, p in enumerate(powers):
                final_coeff *= partial_coeffs[idim][p]
            new_gaussians.append(CartesianGaussian(newcenter.center, newcenter.exp,
                                                   powers=powers, coeff=final_coeff))

        if len(new_gaussians) == 0:  # no non-zero components, just return a zeroed gaussian
            return CartesianGaussian(newcenter, newcenter.exp,
                                     self.powers + other.powers, coeff=0.0)
        elif len(new_gaussians) == 1:
            return new_gaussians[0]
        else:
            return new_gaussians

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

        integ = super().integral
        for p in self.powers:
            if p == 0:  # no contribution
                continue
            elif p % 2 == 1:  # integral of odd function is exactly 0
                return 0.0
            elif p < 0:
                raise ValueError('Powers must be positive or 0')
            else:
                integ *= ODD_FACTORIAL[p-1]/((2.0 * self.exp) ** (p//2))
        return integ
