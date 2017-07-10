""" Tests internal math routines
"""

import random

import itertools
import numpy as np
import numpy.testing as npt
import pytest

import moldesign
from moldesign import units as u

from . import helpers

registered_types = {}


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@typedfixture('gaussian')
def std_1d_gaussian():
    g = moldesign.orbitals.gaussians.Gaussian([0.0]*u.angstrom,
                                              1.0/u.angstrom ** 2)
    return g


@typedfixture('gaussian')
def std_3d_gaussian():
    g = moldesign.orbitals.gaussians.Gaussian([0.0, 0.0, 0.0]*u.angstrom,
                                              1.0/u.angstrom ** 2)
    return g


@typedfixture('gaussian')
def cart_3d_gaussian():
    g = moldesign.orbitals.CartesianGaussian(
            center=[random.random() for i in range(3)]*u.angstrom,
            powers=[1, 3, 0],
            exp=1.1/u.angstrom ** 2,
            coeff=1.0)
    return g


@typedfixture('gaussian')
def spherical_3d_gaussian():
    g = moldesign.orbitals.SphericalGaussian(
            center=[random.random() for i in range(3)]*u.angstrom,
            l=3, m=-2,
            exp=1.1/u.angstrom ** 2,
            coeff=1.0)
    return g


@pytest.mark.parametrize('objkey', registered_types['gaussian'])
@pytest.mark.screening
def test_gaussian_integral_and_dimensionality(objkey, request):
    g = request.getfixturevalue(objkey)
    assert g.ndims == len(g.center)

    intval = g.integral
    expectval = g.coeff*(np.pi/g.exp) ** (g.ndims/2.0)
    _assert_almost_equal(intval,
                         expectval,
                         decimal=10)


def _make_rando_gaussian(powers, withunits=True):
    if withunits:
        length = u.angstrom
    else:
        length = 1.0
    return moldesign.orbitals.CartesianGaussian((np.random.rand(3)-0.5)*1.0 * length,
                                                (random.random()*5)/(length ** 2),
                                                powers=powers,
                                                coeff=random.random())

test_powers = ((0,0,0), (1,0,0), (0,1,0), (0,0,1), (2,0,0), (0,2,0), (0,0,2),
               (1,1,1), (2,2,2), (1,2,3), (0,1,4))
cartesian_test_suite = list(itertools.product(test_powers, test_powers, [True, False]))
cartesian_test_ids = ['[%d%d%d]*[%d%d%d]/%s' % (p[0] + p[1] + ('units' if p[2] else 'c-num',))
                      for p in cartesian_test_suite]

@pytest.mark.parametrize('p1,p2,withunits',
                         cartesian_test_suite,
                         ids=cartesian_test_ids)
def test_gaussian_multiplication_amplitudes(p1, p2, withunits):
    """ Tests that ``g1(x) * g2(x) == (g1 * g2)(x)``
    """
    gauss1 = _make_rando_gaussian(p1, withunits)
    gauss2 = _make_rando_gaussian(p2, withunits)

    testcoords = 6.0*(np.random.rand(50, 3)-0.5)
    if withunits: testcoords = testcoords * u.angstrom

    # test the product of the two AND the squares of each
    #for g1, g2 in itertools.product((gauss1, gauss2), (gauss1, gauss2)):
    g1, g2 = gauss1, gauss2  # temp
    g1g2 = g1 * g2
    if isinstance(g1g2, moldesign.orbitals.AbstractFunction):
        gvals = g1g2(testcoords)
        assert g1g2.coeff != 0.0
    else:
        gvals = sum(gg(testcoords) for gg in g1g2)
    g1vals = g1(testcoords)
    g2vals = g2(testcoords)

    prodvals = g1vals * g2vals

    helpers.assert_almost_equal(prodvals, gvals)


@pytest.mark.parametrize('objkey',
                         registered_types['gaussian'])
def test_gaussian_function_values(objkey, request):
    g = request.getfixturevalue(objkey)

    for idim in range(g.ndims):
        coord = g.center.copy()
        randoffset = 4.0 * (random.random() - 0.5) * g.exp**-0.5
        coord[idim] += randoffset
        funcval = _gfuncval(g, coord)
        retval = g(coord)
        _assert_almost_equal(funcval, retval)


@pytest.mark.parametrize('objkey',
                         registered_types['gaussian'])
@pytest.mark.screening
def test_vectorized_gaussian_function_evaluations(objkey, request):
    g = request.getfixturevalue(objkey)

    coords = np.zeros((5, g.ndims)) * g.center.units
    for i in range(5):
        coords[i] = g.center
        randoffset = 4.0 * (random.random() - 0.5) * g.exp**-0.5
        idim = random.randrange(g.ndims)
        coords[i, idim] += randoffset

    vector_results = g(coords)
    expected = u.array([g(c) for c in coords])
    if vector_results.dimensionless:
        vector_results = vector_results._magnitude

    _assert_almost_equal(vector_results, expected, decimal=8)


def test_normalized_gaussian_self_overlap_is_unity():
    g1 = moldesign.orbitals.gaussians.Gaussian([0.0, 0.0, 0.0]*u.angstrom,
                                              1.0/u.angstrom ** 2,
                                              coeff=-1.0)
    g2 = moldesign.orbitals.gaussians.Gaussian([0.0, 0.0, 0.0]*u.angstrom,
                                              1.0/u.angstrom ** 2,
                                              coeff=123.456)
    npt.assert_almost_equal(-1.0,
                            g1.overlap(g2, normalized=True))


def _gfuncval(g, coord):
    fv = g.coeff * np.exp(-g.exp * np.sum((g.center - coord)**2))
    for r, r0, pow in zip(coord, g.center, g.powers):
        if pow != 0:
            fv *= (r-r0)**pow
    return fv


def _assert_almost_equal(a, b, *args, **kwargs):
    a_is_quantity = hasattr(a,'units')
    b_is_quantity = hasattr(b,'units')

    if not (a_is_quantity or b_is_quantity):
        return np.testing.assert_almost_equal(a, b,
                                              *args, **kwargs)
    else:
        assert a_is_quantity and b_is_quantity
        units = a.units
        return np.testing.assert_almost_equal(a.value_in(units),
                                              b.value_in(units),
                                              *args, **kwargs)


def test_cart_to_powers():
    from moldesign.orbitals.gaussians import cart_to_powers
    assert cart_to_powers('y') == [0, 1, 0]
    assert cart_to_powers('xxyz') == [2, 1, 1]
    assert cart_to_powers('zx^3') == [3,0,1]
    

@pytest.mark.parametrize('key', ['std_3d_gaussian', 'cart_3d_gaussian', 'cart_3d_gaussian'])
def test_numerical_vs_analytical_norm(key, request):
    g = request.getfixturevalue(key)
    ananorm = g.norm

    width = np.sqrt(1/g.exp)
    ranges = np.ones((3,2)) * 5.0 * width
    ranges[:,0] *= -1
    ranges += g.center[:, None]
    grid = moldesign.mathutils.VolumetricGrid(*ranges, 25)
    allpoints = grid.allpoints()

    with np.errstate(under='ignore'):
        vals = g(allpoints)

    numnorm = np.sqrt(grid.dx**3 * (vals**2).sum())
    helpers.assert_almost_equal(ananorm, numnorm)


