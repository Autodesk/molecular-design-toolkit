from __future__ import print_function, absolute_import, division
from future import standard_library
standard_library.install_aliases()
from future.builtins import *

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

__all__ = ['SPHERE_TO_CART']


def cart_to_polar_angles(coords):
    if len(coords.shape) == 2:
        r_xy_2 = np.sum(coords[:,:2]**2, axis=1)
        theta = np.arctan2(np.sqrt(r_xy_2), coords[:,2])
        phi = np.arctan2(coords[:,1], coords[:,0])
        return theta, phi
    else:
        assert len(coords) == 3 and len(coords.shape) == 1
        r_xy_2 = np.sum(coords[:2]**2)
        theta = np.arctan2(np.sqrt(r_xy_2), coords[2])
        phi = np.arctan2(coords[1], coords[0])
        return theta, phi


class Y(object):
    r""" A real-valued spherical harmonic function

    These functions are orthonormalized over the sphere such that

    .. math::
       \int d^2\Omega \: Y_{lm}(\theta, \phi) Y_{l'm'}(\theta, \phi) = \delta_{l,l'} \delta_{m,m'}

    References:
        https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
    """
    def __init__(self, l, m):
        self.l = l
        self.m = m

        self._posneg = -1
        if self.m < 0:
            self._posneg *= -1
        if self.m % 2 == 0:
            self._posneg *= -1

    def __call__(self, coords):
        from scipy.special import sph_harm

        theta, phi = cart_to_polar_angles(coords)
        if self.m == 0:
            return (sph_harm(self.m, self.l, phi, theta)).real

        vplus = sph_harm(abs(self.m), self.l, phi, theta)
        vminus = sph_harm(-abs(self.m), self.l, phi, theta)
        value = np.sqrt(1/2.0) * (self._posneg*vplus + vminus)

        if self.m < 0:
            return -value.imag
        else:
            return value.real


class Cart(object):
    def __init__(self, px, py, pz, coeff):
        self.powers = np.array([px, py, pz], dtype='int')
        self.coeff = coeff
        self._l = self.powers.sum()

    def __call__(self, coords):
        """ Evaluate this function at the given list of coordinates

        Args:
            coords (Vector[len=3] or Matrix[matrix=shape(*,3)): coordinate(s) at which to evaluate

        Returns:
            Scalar or Vector: value of the function at these coordinates
        """

        c = coords**self.powers

        if len(coords.shape) == 2:
            r_l = (coords**2).sum(axis=1) ** (-self._l / 2.0)
            return self.coeff * np.product(c, axis=1) * r_l
        else:
            r_l = (coords**2).sum() ** (-self._l / 2.0)
            return self.coeff * np.product(c) * r_l

    def __iter__(self):
        yield self  # for interface compatibility with the CartSum class


class CartSum(object):
    def __init__(self, coeff, carts):
        self.coeff = coeff

        prefacs = []
        powers = []
        self.l = None
        for cart in carts:
            powers.append(np.array(cart[:3], dtype='int'))
            prefacs.append(np.array(cart[3]))
            if self.l is None:
                self.l = sum(cart[:3])
            else:
                assert self.l == sum(cart[:3])

        self.prefactors = np.array(prefacs)
        self.powers = np.array(powers)

    def __call__(self, coords):
        # likely can speed this up a lot using clever numpy broadcasting

        if len(coords.shape) == 2:
            c = np.zeros((len(coords), ))
            axis = 1
            r_l = (coords**2).sum(axis=1) ** (-self.l / 2.0)
        else:
            c = 0.0
            axis = None
            r_l = (coords**2).sum() ** (-self.l / 2.0)

        for factor, power in zip(self.prefactors, self.powers):
            c += factor * np.product(coords**power, axis=axis)

        return c * self.coeff * r_l

    def __iter__(self):
        for pf, (px, py, pz) in zip(self.prefactors, self.powers):
            yield Cart(px, py, pz, coeff=self.coeff*pf)



def sqrt_x_over_pi(num, denom):
    return np.sqrt(num / (denom*np.pi))


                  ############   s     ############
SPHERE_TO_CART = {(0, 0): Cart(0, 0, 0, sqrt_x_over_pi(1, 4)),

                  ############   p     ############
                  (1, -1): Cart(0, 1, 0, sqrt_x_over_pi(3, 4)),
                  (1, 0): Cart(0, 0, 1, sqrt_x_over_pi(3, 4)),
                  (1, 1): Cart(1, 0, 0, sqrt_x_over_pi(3, 4)),

                  ############   d     ############
                  (2, -2): Cart(1, 1, 0, sqrt_x_over_pi(15, 4)),
                  (2, -1): Cart(0, 1, 1, sqrt_x_over_pi(15, 4)),
                  (2, 0): CartSum(sqrt_x_over_pi(5, 16),
                                  [(2, 0, 0, -1.0),
                                   (0, 2, 0, -1.0),
                                   (0, 0, 2, 2.0)]),
                  (2, 1): Cart(1, 0, 1, sqrt_x_over_pi(15, 4)),
                  (2, 2): CartSum(sqrt_x_over_pi(15, 16),
                                  [(2, 0, 0, 1.0),
                                   (0, 2, 0, -1.0)]),

                  ############   f     ############
                  (3, -3): CartSum(sqrt_x_over_pi(35, 32),
                                   [(2, 1, 0, 3.0),
                                    (0, 3, 0, -1.0)]),

                  (3, -2): Cart(1, 1, 1, sqrt_x_over_pi(105, 4)),

                  (3, -1): CartSum(sqrt_x_over_pi(21, 32),
                                   [(0, 1, 2, 4.0),
                                    (2, 1, 0, -1.0),
                                    (0, 3, 0, -1.0)]),

                  (3, 0): CartSum(sqrt_x_over_pi(7, 16),
                                  [(0, 0, 3, 2.0),
                                   (2, 0, 1, -3.0),
                                   (0, 2, 1, -3.0)]),

                  (3, 1): CartSum(sqrt_x_over_pi(21, 32),
                                   [(1, 0, 2, 4.0),
                                    (3, 0, 0, -1.0),
                                    (1, 2, 0, -1.0)]),

                  (3, 2): CartSum(sqrt_x_over_pi(105, 16),
                                  [(2, 0, 1, 1.0),
                                   (0, 2, 1, -1.0)]),

                  (3, 3): CartSum(sqrt_x_over_pi(35, 32),
                                  [(3, 0, 0, 1.0),
                                   (1, 2, 0, -3.0)]),
                  }