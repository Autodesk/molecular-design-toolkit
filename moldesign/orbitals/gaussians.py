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


class AbstractFunction(object):
    """ Abstract base class for basis functions
    """
    def normalize(self):
        """ Give this function unit norm by adjusting its coefficient
        """
        self.coeff /= np.sqrt(self.norm)

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
            integral /= np.sqrt(self.norm*other.norm)
        return integral

    @property
    def norm(self):
        r""" The L2-Norm of this gaussian:

        .. math::
            \sqrt{\int \left| G(\mathbf r) \right|^2 d^N \mathbf r}
        """
        return np.sqrt(self.overlap(self))

    def __str__(self):
        return "%d-D %s with norm %s" % (self.ndims, self.__class__, self.norm)


class CartesianGaussian(AbstractFunction):
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
        self.center = center
        self.exp = exp
        self.powers = np.array(powers)
        if coeff is None:
            self.coeff = 1.0  # dummy value overwritten by self.normalize()
            self.normalize()
        else:
            self.coeff = coeff

        self.shell = sum(self.powers)

    def __repr__(self):
        return ("<{ndim}-D Gaussian (cart) (coeff: {coeff:4.2f}, "
                "cartesian powers: {powers}, "
                "exponent: {exp:4.2f}, "
                "center: {center}>").format(
                ndim=self.ndim,
                center=self.center, exp=self.exp,
                powers=tuple(self.powers), coeff=self.coeff)

    @property
    def ndim(self):
        return len(self.powers)
    num_dimensions = ndims = ndim

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
            Scalar: function value(s) at the passed coordinates
        """
        if len(coords.shape) > 1:
            axis = 1
        else:
            axis = None

        disp = coords - self.center
        prd = disp*disp  # don't use np.dot - allow multiple coords at once
        r2 = prd.sum(axis=axis)

        result = self.coeff * np.exp(-self.exp * r2)
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
        return np.product(disp.magnitude**self.powers, axis=axis)*disp.units ** self.powers.sum()

    def __mul__(self, other):
        """ Returns product of two gaussian functions, which is also a gaussian

        Args:
            other (CartesianGaussian): other gaussian wavepacket

        Returns:
            CartesianGaussian: product gaussian
        """
        if (self.center != other.center).any():
            raise NotImplementedError()

        # convert widths to prefactor form
        a = self.exp
        b = other.exp
        exp = a + b
        center = (a*self.center + b*other.center)/(a+b)
        powers = self.powers + other.powers
        return CartesianGaussian(center=center, exp=exp,
                                 powers=powers, coeff=self.coeff*other.coeff)

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
        integ = (np.pi/self.exp)**(self.ndim/2.0)
        for p in self.powers:
            if p == 0:  # no contribution
                continue
            elif p % 2 == 1:  # integral of odd function is exactly 0
                return 0.0
            elif p < 0:
                raise ValueError('Powers must be positive or 0')
            else:
                integ *= _ODD_FACTORIAL[p-1]/(2 ** (p+1))
        return self.coeff * integ


def Gaussian(center, exp, coeff=1.0):
    r""" Constructor for an N-dimensional gaussian function.

    The function is given by:
    .. math::
        G(\mathbf r) = C e^{-a\left| \mathbf r - \mathbf r_0 \right|^2}

    where *C* is ``self.coeff``, *a* is ``self.exp``, and :math:`\mathbf r_0` is ``self.center``.

    Note:
        This is just a special case of a cartesian gaussian where all the powers are 0.
    """
    return CartesianGaussian(center, exp,
                             powers=[0 for x in center],
                             coeff=coeff)


class SphericalGaussian(AbstractFunction):
    r""" Stores a 3-dimensional real spherical gaussian function:

    .. math::
        G_{nlm}(\mathbf r) = C Y^l_m(\mathbf r - \mathbf r_0) r^l e^{-a\left| \mathbf r - \mathbf r_0 \right|^2}

    where *C* is ``self.coeff``, *a* is ``self.exp``, and :math:`\mathbf r_0` is ``self.center``,
    *(l,m)* are given by ``(self.l, self.m)``, and :math:`Y^l_m(\mathbf{r})` are the _real_
    spherical harmonics.

    References:
        Schlegel and Frisch. Transformation between cartesian and pure spherical harmonic gaussians.
            Int J Quantum Chem 54, 83-87 (1995). doi:10.1002/qua.560540202

        https://en.wikipedia.org/wiki/Table_of_spherical_harmonics#Real_spherical_harmonics
    """

    ndims = 3

    def __init__(self, center, exp, l, m, coeff=None):
        from moldesign.mathutils import spherical_harmonics

        self.center = center
        assert len(self.center) == 3
        self.exp = exp
        self.l = l
        self.m = m

        if coeff is None:
            self.coeff = 1.0
            self.normalize()
        else:
            self.coeff = coeff

        self._spherical_harmonic = spherical_harmonics.Y(l, m)

    def __repr__(self):
        return ("<3D Gaussian (Spherical) (coeff: {coeff:4.2f}, "
                "exponent: {exp:4.2f}, "
                "(l,m) = {qnums}").format(
                center=self.center, exp=self.exp, coeff=self.coeff,
                qnums=(self.l, self.m))

    def __call__(self, coords, _include_angular=True):
        """ Evaluate this function at the given coordinates.

        Can be called either with a 1D column (e.g., ``[1,2,3]*u.angstrom ``) or
        an ARRAY of coordinates (``[[0,0,0],[1,1,1]] * u.angstrom``)

        Args:
            _include_angular (bool): include the contribution from the non-exponential parts
                (for computational efficiency, this can be omitted now and included later)

        Args:
            coords (Vector[length]): 3D Coordinates or list of 3D coordinates

        Returns:
            Scalar: function value(s) at the passed coordinates
        """
        if len(coords.shape) > 1:
            axis = 1
        else:
            axis = None

        disp = coords - self.center
        prd = disp*disp
        r2 = prd.sum(axis=axis)

        result = self.coeff * np.exp(-self.exp * r2)
        if _include_angular:
            result *= r2**(self.l/2.0) * self._spherical_harmonic(disp)
        return result

    def angular_part(self, coords):
        if len(coords.shape) > 1:
            axis = 1
        else:
            axis = None

        disp = coords - self.center
        prd = disp*disp
        r2 = prd.sum(axis=axis)

        return r2**(self.l/2.0) * self._spherical_harmonic(disp)


# Precompute odd factorial values (N!!)
_ODD_FACTORIAL = {0: 1}  # by convention
_ofact = 1
for _i in range(1, 20, 2):
    _ofact *= _i
    _ODD_FACTORIAL[_i] = float(_ofact)


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
