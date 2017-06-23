import pytest

import numpy as np

from moldesign import mathutils
from moldesign import units as u


def y_numpy_unit_vector():
    v = np.zeros(3)
    v[1] = 1.0
    return v


def y_length_vector():
    v = y_numpy_unit_vector() * 3.0 * u.nm
    return v


@pytest.mark.parametrize('y_unit_vector', [y_numpy_unit_vector, y_length_vector])
def test_perpendicular_to_y(y_unit_vector):
    arr = y_unit_vector()

    outvec = mathutils.perpendicular(arr)
    assert isinstance(outvec, np.ndarray)

    assert outvec[1] == 0.0
    assert abs(mathutils.norm(outvec) - 1.0) < 1e-12


def test_norm_numpy_y_unit_vector():
    arr = y_numpy_unit_vector()
    assert mathutils.norm(arr) == 1.0
    assert (mathutils.normalized(arr) == arr).all()


def test_normalization_with_units():
    arr = y_length_vector()
    norm = mathutils.norm(arr)
    assert arr.units == norm.units
    assert (arr/norm).dimensionless

    normalized = mathutils.normalized(arr)
    assert isinstance(normalized, np.ndarray) or normalized.units.dimensionless
    np.testing.assert_allclose(normalized,
                               arr/norm)
    assert abs(mathutils.norm(normalized) - 1.0) < 1e-12


