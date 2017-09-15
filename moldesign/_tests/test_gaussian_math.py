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
from .molecule_fixtures import *

registered_types = {}

__PYTEST_MARK__ = ['math', 'gaussians']


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""
    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@pytest.fixture
def std_1d_gaussian():
    g = moldesign.orbitals.gaussians.Gaussian([0.0]*u.angstrom,
                                              1.0/u.angstrom ** 2)
    return g


@typedfixture('basis_fn')
def std_3d_gaussian():
    g = moldesign.orbitals.gaussians.Gaussian([0.0, 0.0, 0.0]*u.angstrom,
                                              1.0/u.angstrom ** 2)
    return g


@typedfixture('basis_fn')
def cartesian_3d_gaussian():
    g = moldesign.orbitals.CartesianGaussian(
            center=[random.random() for i in range(3)]*u.angstrom,
            powers=[1, 3, 0],
            alpha=1.1/u.angstrom ** 2,
            coeff=1.0)
    return g


@typedfixture('basis_fn')
def spherical_3d_gaussian():
    g = moldesign.orbitals.SphericalGaussian(
            center=[random.random() for i in range(3)]*u.angstrom,
            l=3, m=-2,
            alpha=1.1/u.angstrom ** 2,
            coeff=1.0)
    return g


@pytest.mark.parametrize('objkey', ['std_1d_gaussian','std_3d_gaussian'])
@pytest.mark.screening
def test_gaussian_integral_and_dimensionality(objkey, request):
    g = request.getfixturevalue(objkey)
    assert g.ndim == len(g.center)

    intval = g.integral
    expectval = g.coeff*(np.pi/g.alpha) ** (g.ndim/2.0)
    _assert_almost_equal(intval,
                         expectval,
                         decimal=10)


@pytest.fixture
def linear_combination():
    return _make_rando_linear_combination(True)


def _make_rando_gaussian(withunits=True):
    if withunits:
        length = u.angstrom
    else:
        length = 1.0
    return moldesign.orbitals.Gaussian((np.random.rand(3)-0.5)*1.0 * length,
                                       (random.random()*5)/(length ** 2),
                                       coeff=random.random())


def _make_rando_cartesian_gaussian(powers, withunits=True):
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


def _make_rando_linear_combination(withunits=True):
    gaussians = []
    if withunits:
        length = u.angstrom
    else:
        length = 1.0

    center = (np.random.rand(3)-0.5)*1.0 * length

    for pwr in [(0,0,0), (1,1,1), (3,2,1)]:
        gaussians.append(
                moldesign.orbitals.CartesianGaussian(
                        center=center,
                        powers=pwr,
                        alpha=(10.0 * (random.random()+3))/(length**2),
                        coeff=1/(np.sqrt(3.0))))
    lc = moldesign.orbitals.PrimitiveSum(gaussians)
    lc.ndims = 3  # so it works with the test suite
    return lc


@pytest.mark.parametrize('withunits', [True, False], ids=['quantity','number'])
def test_numerical_vs_analytical_overlap_gauss(withunits):
    p1 = _make_rando_gaussian(withunits)
    p2 = _make_rando_gaussian(withunits)
    _assert_numerical_analytical_overlaps_match(p1, p2)


@pytest.mark.parametrize('withunits', [True, False], ids=['quantity','number'])
def test_numerical_vs_analytical_overlap_cartesian(withunits):
    p1 = _make_rando_cartesian_gaussian((1,2,3), withunits)
    p2 = _make_rando_cartesian_gaussian((1,0,1), withunits)
    _assert_numerical_analytical_overlaps_match(p1, p2)


@pytest.mark.parametrize('withunits', [True, False], ids=['quantity','number'])
def test_numerical_vs_analytical_overlap_spherical(withunits):
    p1 = _make_rando_spherical_gaussian(1,-1, withunits)
    p2 = _make_rando_spherical_gaussian(2,0, withunits)
    _assert_numerical_analytical_overlaps_match(p1, p2)


@pytest.mark.parametrize('withunits', [True, False], ids=['quantity','number'])
def test_numerical_vs_analytical_overlap_linear_combination(withunits):
    p1 = _make_rando_linear_combination(withunits)
    p2 = _make_rando_linear_combination(withunits)
    _assert_numerical_analytical_overlaps_match(p1, p2)


def _assert_numerical_analytical_overlaps_match(g1, g2):
    olap = g1.overlap(g2)
    try:
        prod = g1*g2
    except NotImplementedError:
        assert isinstance(g1, moldesign.orbitals.SphericalGaussian)
        assert isinstance(g2, moldesign.orbitals.SphericalGaussian)
    else:
        helpers.assert_almost_equal(prod.integral, olap)

    def assert_with_resolution(npoints):
        allpoints, grid = helpers.generate_grid(g1, g2, npoints)
        with np.errstate(under='ignore'):
            prodvals = g1(allpoints) * g2(allpoints)
        numsum = prodvals.sum() * grid.dx * grid.dy * grid.dz
        helpers.assert_almost_equal(numsum, olap, decimal=4)

    # If numerical isn't equal to analytical, try again with higher resolution
    # to make sure the failure isn't due to a sparse grid:
    try:
        assert_with_resolution(64)
    except AssertionError:
        pass
    else:
        return

    try:
        assert_with_resolution(128)
    except AssertionError:
        pass
    else:
        return

    assert_with_resolution(256)




@pytest.mark.parametrize('withunits', [False, True])
def test_gaussian_multiplication_amplitudes(withunits):
    g1 = _make_rando_gaussian(withunits)
    g2 = _make_rando_gaussian(withunits)
    _assert_same_function_values(g1, g2, withunits)


# parameterizations across a sample of cartesian gaussians
test_powers = ((0,0,0), (1,0,0), (0,1,0), (0,0,1), (2,0,0), (1,1,1), (2,0,2), (4,1,1))
cartesian_test_suite = list(itertools.product(test_powers, test_powers, [True, False]))
cartesian_test_ids = ['[%d%d%d]*[%d%d%d]/%s' % (p[0] + p[1] + ('units' if p[2] else 'c-num',))
                      for p in cartesian_test_suite]


@pytest.mark.parametrize('p1,p2,withunits',
                         cartesian_test_suite,
                         ids=cartesian_test_ids)
def test_cartesian_gaussian_multiplication_amplitudes(p1, p2, withunits):
    """ Tests that ``g1(x) * g2(x) == (g1 * g2)(x)``
    """
    g1 = _make_rando_cartesian_gaussian(p1, withunits)
    g2 = _make_rando_cartesian_gaussian(p2, withunits)
    _assert_same_function_values(g1, g2, withunits)


def _assert_same_function_values(g1, g2, withunits):
    testcoords = 6.0*(np.random.rand(50, 3)-0.5)
    if withunits:
        testcoords = testcoords*u.angstrom
    g1g2 = g1*g2
    gvals = g1g2(testcoords)
    g1vals = g1(testcoords)
    g2vals = g2(testcoords)
    prodvals = g1vals*g2vals
    helpers.assert_almost_equal(prodvals, gvals)


def test_initial_gaussian_normalization_gaussian():
    center = np.random.rand(3) * u.angstrom
    exp = 5.12 / u.angstrom**2
    g2 = moldesign.orbitals.Gaussian(center, exp, normalized=True)
    helpers.assert_almost_equal(1.0, _numerical_norm(g2), decimal=3)
    helpers.assert_almost_equal(1.0, g2.norm)


def test_initial_gaussian_normalization_with_prefactor():
    center = np.random.rand(3) * u.angstrom
    exp = 5.12 / u.angstrom**2
    g1 = moldesign.orbitals.Gaussian(center, exp, coeff=3.0*u.angstrom, normalized=True)
    helpers.assert_almost_equal(3.0*u.angstrom, _numerical_norm(g1), decimal=3)
    helpers.assert_almost_equal(3.0*u.angstrom, g1.norm)


def test_initial_normalization_cartesian():
    center = np.random.rand(3) * u.angstrom
    exp = 5.12 / u.angstrom**2
    for powers in itertools.product(range(4), range(4), range(4)):
        g2 = moldesign.orbitals.CartesianGaussian(center, exp, powers, normalized=True)
        helpers.assert_almost_equal(1.0, _numerical_norm(g2), decimal=3)
        helpers.assert_almost_equal(1.0, g2.norm)


def test_initial_normalization_cartesian_with_prefactor():
    center = np.random.rand(3) * u.angstrom
    exp = 5.12 / u.angstrom**2
    for powers in itertools.product(range(4), range(4), range(4)):
        g1 = moldesign.orbitals.CartesianGaussian(center, exp, powers, coeff=3.0, normalized=True)
        helpers.assert_almost_equal(3.0, _numerical_norm(g1), decimal=3)
        helpers.assert_almost_equal(3.0, g1.norm)


def test_initial_normalization_spherical():
    center = np.random.rand(3) * u.angstrom
    exp = 5.12 / u.angstrom**2
    for l in range(5):
        for m in range(-l, l+1):
            g2 = moldesign.orbitals.SphericalGaussian(center, exp, l, m, normalized=True)
            helpers.assert_almost_equal(1.0, _numerical_norm(g2), decimal=3)
            helpers.assert_almost_equal(1.0, g2.norm)


def test_initial_normalization_spherical_with_prefactor():
    center = np.random.rand(3) * u.angstrom
    exp = 5.12 / u.angstrom**2
    for l in range(5):
        for m in range(-l, l+1):
            g1 = moldesign.orbitals.SphericalGaussian(center, exp, l, m,
                                                      coeff=3.0 * u.angstrom, normalized=True)
            helpers.assert_almost_equal(3.0 * u.angstrom, _numerical_norm(g1), decimal=3)
            helpers.assert_almost_equal(3.0 * u.angstrom, g1.norm)


def _numerical_norm(g):
    allpoints, grid = helpers.generate_grid(g)
    with np.errstate(under='ignore'):
        vals = g(allpoints)
    numnorm = np.sqrt(grid.dx * grid.dy * grid.dz * (vals**2).sum())
    return numnorm


@pytest.mark.parametrize('objkey', registered_types['basis_fn'])
def test_gaussian_function_values(objkey, request):
    g = request.getfixturevalue(objkey)

    for idim in range(g.ndims):
        coord = g.center.copy()
        randoffset = 4.0 * (random.random() - 0.5) * g.alpha**-0.5
        coord[idim] += randoffset
        funcval = _gfuncval(g, coord)
        retval = g(coord)
        _assert_almost_equal(funcval, retval)


@pytest.mark.parametrize('objkey', registered_types['basis_fn'])
def test_vectorized_gaussian_function_evaluations(objkey, request):
    g = request.getfixturevalue(objkey)

    coords = np.zeros((5, g.ndims)) * g.center.units
    for i in range(5):
        coords[i] = g.center
        randoffset = 4.0 * (random.random() - 0.5) * g.alpha**-0.5
        idim = random.randrange(g.ndims)
        coords[i, idim] += randoffset

    vector_results = g(coords)
    expected = u.array([g(c) for c in coords])
    if vector_results.dimensionless:
        vector_results = vector_results._magnitude

    _assert_almost_equal(vector_results, expected, decimal=8)


@pytest.mark.parametrize('objkey', registered_types['basis_fn'] + ['linear_combination'])
def test_gaussian_str_and_repr_works(objkey, request):
    g1 = request.getfixturevalue(objkey)
    str(g1)
    repr(g1)


@pytest.mark.parametrize('objkey', registered_types['basis_fn'])
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


@pytest.mark.parametrize('objkey', registered_types['basis_fn'])
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


def test_linear_combination_normalization(linear_combination):
    g1 = linear_combination
    oldnorm = g1.norm

    prefactor = (random.random() - 0.5) * 428.23
    for prim in g1:
        prim.coeff *= prefactor

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
    fv = g.coeff * np.exp(-g.alpha * r2)
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
    

@pytest.mark.parametrize('key', registered_types['basis_fn'] + ['linear_combination'])
def test_numerical_vs_analytical_norm(key, request):
    g = request.getfixturevalue(key)
    numnorm = _numerical_norm(g)
    helpers.assert_almost_equal(g.norm, numnorm)


@pytest.mark.screening
def test_s_orbitals_equivalent_among_gaussian_types():
    center = np.random.rand(3) * u.angstrom
    exp = 5.12 / u.angstrom**2

    g_bare = moldesign.orbitals.Gaussian(center, exp)
    g_cart = moldesign.orbitals.CartesianGaussian(center, exp, [0,0,0])
    g_sphr = moldesign.orbitals.SphericalGaussian(center, exp, 0, 0)

    for gauss in (g_bare, g_cart, g_sphr):
        # normalize to amplitude of 1.0 at center
        gauss.coeff = gauss.coeff / gauss(center)
        assert gauss(center) == 1.0

    _assert_orbitals_equivalent(g_bare, g_cart)
    _assert_orbitals_equivalent(g_bare, g_sphr)


# real spherical harmonics that can be represented as a single cartesian term:
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

    g_cart = moldesign.orbitals.CartesianGaussian(center, exp, powers)
    g_sphr = moldesign.orbitals.SphericalGaussian(center, exp, *lm)

    _assert_orbitals_equivalent(g_cart, g_sphr)


@pytest.mark.parametrize('l', range(4), ids=lambda x:'l=%s' % x)
def test_spherical_to_cartesian(l):
    for m in range(-l,l+1):
        center = np.random.rand(3)*u.angstrom
        exp = random.random()*2.0/u.angstrom ** 2

        bf = moldesign.orbitals.SphericalGaussian(center, exp, l, m, normalized=True)
        _assert_orbitals_equivalent(bf, bf.to_cart())


def _assert_orbitals_equivalent(g1, g2):
    helpers.assert_almost_equal(g1.norm,
                                g2.norm)

    testcoords = 6.0*(np.random.rand(50, 3)-0.5)*u.angstrom
    helpers.assert_almost_equal(g1(testcoords),
                                g2(testcoords))


def test_pyscf_and_mdt_norms_are_the_same(h2_rhf_augccpvdz):
    mol = h2_rhf_augccpvdz
    basis = mol.wfn.aobasis

    for bf in basis:
        assert abs(bf.norm - 1.0) < 1e-12


def test_pyscf_and_mdt_overlaps_are_the_same(h2_rhf_augccpvdz):
    mol = h2_rhf_augccpvdz
    basis = mol.wfn.aobasis

    calc_overlap_mat = []
    for i in range(len(basis)):
        calc_overlap_mat.append(
                [basis[i].overlap(basis[j]) for j in range(len(basis))]
        )

    overlaps = u.array(calc_overlap_mat)
    assert isinstance(overlaps, np.ndarray) or overlaps.units == u.dimensionless

    np.testing.assert_allclose(mol.wfn.aobasis.overlaps,
                               overlaps,
                               atol=5.0e-7)
