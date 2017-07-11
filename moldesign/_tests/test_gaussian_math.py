""" Tests internal math routines
"""

import random

import itertools
import numpy as np
import pytest

import moldesign
from moldesign import units as u
from moldesign.mathutils import spherical_harmonics as harmonics

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


@typedfixture('cart-gaussian')
def cart_3d_gaussian():
    g = moldesign.orbitals.CartesianGaussian(
            center=[random.random() for i in range(3)]*u.angstrom,
            powers=[1, 3, 0],
            exp=1.1/u.angstrom ** 2,
            coeff=1.0)
    return g


@typedfixture('spherical-gaussian')
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


def _make_rando_cart_gaussian(powers, withunits=True):
    if withunits:
        length = u.angstrom
    else:
        length = 1.0
    return moldesign.orbitals.CartesianGaussian((np.random.rand(3)-0.5)*1.0 * length,
                                                (random.random()*5)/(length ** 2),
                                                powers=powers,
                                                coeff=random.random())


def _make_rando_spherical_gaussian(l,m, withunits=True):
    if withunits:
        length = u.angstrom
    else:
        length = 1.0
    return moldesign.orbitals.SphericalGaussian((np.random.rand(3)-0.5)*1.0 * length,
                                                (random.random()*5)/(length ** 2),
                                                l,m,
                                                coeff=random.random())

# parameterizations across a sample of cartesian gaussians
test_powers = ((0,0,0), (1,0,0), (0,1,0), (0,0,1), (2,0,0), (1,1,1), (2,0,2), (4,1,1))
cartesian_test_suite = list(itertools.product(test_powers, test_powers, [True, False]))
cartesian_test_ids = ['[%d%d%d]*[%d%d%d]/%s' % (p[0] + p[1] + ('units' if p[2] else 'c-num',))
                      for p in cartesian_test_suite]


@pytest.mark.parametrize('p1,p2,withunits',
                         cartesian_test_suite,
                         ids=cartesian_test_ids)
def test_cart_gaussian_multiplication_amplitudes(p1, p2, withunits):
    """ Tests that ``g1(x) * g2(x) == (g1 * g2)(x)``
    """
    g1 = _make_rando_cart_gaussian(p1, withunits)
    g2 = _make_rando_cart_gaussian(p2, withunits)

    testcoords = 6.0*(np.random.rand(50, 3)-0.5)
    if withunits:
        testcoords = testcoords*u.angstrom
    g1g2 = g1*g2
    gvals = g1g2(testcoords)
    g1vals = g1(testcoords)
    g2vals = g2(testcoords)
    prodvals = g1vals*g2vals
    helpers.assert_almost_equal(prodvals, gvals)


@pytest.mark.parametrize('objkey',
                         registered_types['gaussian'] +
                         registered_types['cart-gaussian'] +
                         registered_types['spherical-gaussian'])
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


@pytest.mark.parametrize('objkey',
                         registered_types['gaussian'] +
                         registered_types['cart-gaussian'] +
                         registered_types['spherical-gaussian'])
def test_gaussian_str_and_repr_works(objkey, request):
    g1 = request.getfixturevalue(objkey)
    str(g1)
    repr(g1)


@pytest.mark.parametrize('objkey',
                         registered_types['gaussian'] +
                         registered_types['cart-gaussian'] +
                         registered_types['spherical-gaussian'])
def test_normalized_gaussian_self_overlap_is_unity(objkey, request):
    g1 = request.getfixturevalue(objkey)
    g2 = g1.copy()
    g1.coeff = -10.0
    g2.coeff = 12341.1832
    olap = g1.overlap(g2, normalized=True)
    assert abs(-1.0 - olap) < 1e-12

    g1.coeff = 10.0
    olap = g1.overlap(g2, normalized=True)
    assert abs(1.0 - olap) < 1e-12


@pytest.mark.parametrize('objkey',
                         registered_types['gaussian'] +
                         registered_types['cart-gaussian'] +
                         registered_types['spherical-gaussian'])
def test_normalization(objkey, request):
    g1 = request.getfixturevalue(objkey)
    oldnorm = g1.norm

    g1.coeff = (random.random() - 0.5) * 428.23
    try:
        assert g1.norm != oldnorm
    except u.DimensionalityError:
        pass  # this is a reasonable thing to happen too

    g1.normalize()
    assert abs(g1.norm - 1.0) < 1e-12


def _gfuncval(g, coord):
    r = g.center - coord
    if len(coord.shape) > 1:
        r2 = np.sum(r**2, axis=1)
    else:
        r2 = np.sum(r**2)
    fv = g.coeff * np.exp(-g.exp * r2)
    if isinstance(g, moldesign.orbitals.SphericalGaussian):
        fv *= r2**(g.l/2.0) * harmonics.Y(g.l, g.m)(coord - g.center)
    elif isinstance(g, moldesign.orbitals.CartesianGaussian):  # assume cartesian
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


def test_convert_cartesian_label_to_array_of_integer_powers():
    from moldesign.orbitals.gaussians import cart_to_powers
    assert cart_to_powers('y') == [0, 1, 0]
    assert cart_to_powers('xxyz') == [2, 1, 1]
    assert cart_to_powers('zx^3') == [3,0,1]
    

@pytest.mark.parametrize('key', ['std_3d_gaussian', 'cart_3d_gaussian', 'spherical_3d_gaussian'])
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


@pytest.mark.screening
def test_s_orbitals_equivalent_among_gaussian_types():
    center = np.random.rand(3) * u.angstrom
    exp = 5.12 / u.angstrom**2
    testcoords = 6.0*(np.random.rand(50, 3)-0.5) * u.angstrom

    g_bare = moldesign.orbitals.Gaussian(center, exp)
    g_cart = moldesign.orbitals.CartesianGaussian(center, exp, [0,0,0])
    g_sphr = moldesign.orbitals.SphericalGaussian(center, exp, 0, 0)

    for gauss in (g_bare, g_cart, g_sphr):
        # normalize to amplitude of 1.0 at center
        gauss.coeff = gauss.coeff / gauss(center)
        assert gauss(center) == 1.0

    barevals = g_bare(testcoords)
    cartvals = g_cart(testcoords)
    spherevals = g_sphr(testcoords)

    helpers.assert_almost_equal(barevals, cartvals)
    helpers.assert_almost_equal(barevals, spherevals)

    helpers.assert_almost_equal(g_bare.norm, g_cart.norm)
    helpers.assert_almost_equal(g_sphr.norm, g_cart.norm)


LM_TO_CART = {(1,-1): (0,1,0),
              (1,0): (0,0,1),
              (1,1): (1,0,0),
              (2,-2): (1,1,0),
              (2,-1): (0,1,1),
              (2,1): (1,0,1),
              (3,-2): (1,1,1)}


@pytest.mark.parametrize('lm,powers',
                         LM_TO_CART.items(),
                         ids=['lm:%d,%d, xyz:%d%d%d' % (args[0] + args[1])
                              for args in LM_TO_CART.items()])
def test_orbitals_same_in_cartesian_and_spherical(lm, powers):
    center = np.random.rand(3) * u.angstrom
    exp = 5.12 / u.angstrom**2
    testcoords = 6.0*(np.random.rand(50, 3)-0.5) * u.angstrom

    g_cart = moldesign.orbitals.CartesianGaussian(center, exp, powers)
    g_sphr = moldesign.orbitals.SphericalGaussian(center, exp, *lm)

    helpers.assert_almost_equal(g_cart(testcoords),
                                g_sphr(testcoords))

    helpers.assert_almost_equal(g_sphr.norm, g_cart.norm)
