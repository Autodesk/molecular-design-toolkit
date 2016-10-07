""" Tests internal math routines
"""

import random

import numpy as np
import numpy.testing as npt
import pytest

import moldesign
from moldesign import units as u

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
    g = moldesign.methods.gaussians.Gaussian([0.0]*u.angstrom,
                                             1.0/u.angstrom ** 2)
    return g


@typedfixture('gaussian')
def std_3d_gaussian():
    g = moldesign.methods.gaussians.Gaussian([0.0, 0.0, 0.0]*u.angstrom,
                                             1.0/u.angstrom ** 2)
    return g


@typedfixture('gaussian')
def std_10d_gaussian():
    g = moldesign.methods.gaussians.Gaussian([0.0 for i in xrange(10)]*u.angstrom,
                                             1.0/u.angstrom ** 2)
    return g

@typedfixture('cartesiangaussian')
def cart_10d_gaussian():
    g = moldesign.methods.gaussians.CartesianGaussian(center=[random.random() for i in xrange(10)]*u.angstrom,
                                                      powers=[random.choice([0, 1, 2, 3]) for i in xrange(10)],
                                                      exp=1.1/u.angstrom ** 2)
    return g



@pytest.mark.parametrize('objkey', registered_types['gaussian'])
def test_gaussian_integral_and_dimensionality(objkey, request):
    g = request.getfuncargvalue(objkey)
    assert g.ndims == len(g.center)

    intval = g.integral
    expectval = g.coeff*(np.pi/g.exp) ** (g.ndims/2.0)
    _assert_almost_equal(intval,
                         expectval,
                         decimal=10)


@pytest.mark.parametrize('objkey',
                         registered_types['gaussian'] + registered_types['cartesiangaussian'])
def test_gaussian_function_values(objkey, request):
    g = request.getfuncargvalue(objkey)

    for idim in xrange(g.ndims):
        coord = g.center.copy()
        randoffset = 4.0 * (random.random() - 0.5) * g.exp**-0.5
        coord[idim] += randoffset
        funcval = _gfuncval(g, coord)
        retval = g(coord)
        _assert_almost_equal(funcval, retval)


@pytest.mark.parametrize('objkey',
                         registered_types['gaussian'] + registered_types['cartesiangaussian'])
def test_vectorized_gaussian_function_evaluations(objkey, request):
    g = request.getfuncargvalue(objkey)

    coords = np.zeros((5, g.ndims)) * g.center.units
    for i in xrange(5):
        coords[i] = g.center
        randoffset = 4.0 * (random.random() - 0.5) * g.exp**-0.5
        idim = random.randrange(g.ndims)
        coords[i, idim] += randoffset

    vector_results = g(coords)
    expected = u.array([g(c) for c in coords])

    _assert_almost_equal(vector_results, expected, decimal=8)


def test_gaussian_self_overlap_is_unity():
    g1 = moldesign.methods.gaussians.Gaussian([0.0, 0.0, 0.0]*u.angstrom,
                                              1.0/u.angstrom ** 2,
                                              coeff=-1.0)
    g2 = moldesign.methods.gaussians.Gaussian([0.0, 0.0, 0.0]*u.angstrom,
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
