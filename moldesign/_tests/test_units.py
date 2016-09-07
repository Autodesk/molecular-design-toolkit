""" Tests the unit system
"""
import pytest
import numpy as np

from moldesign import units


def test_scalar_comparison_dimensionality_errors():
    with pytest.raises(units.DimensionalityError):
        x = 1.0 * units.angstrom == 1.0*units.ureg.kilograms

    with pytest.raises(units.DimensionalityError):
        y = 1.0 * units.fs < 1.0 * units.ureg.meter

    with pytest.raises(units.DimensionalityError):
        z = 1.0 * units.ureg.hectare >= 1.0 * units.ureg.astronomical_unit


def test_array_comparison_dimensionality_errors():
    mylist = [0.0, -1.0, 1.0]

    with pytest.raises(units.DimensionalityError):
        x = mylist * units.angstrom == mylist*units.ureg.kilograms

    with pytest.raises(units.DimensionalityError):
        y = mylist * units.fs < mylist * units.ureg.meter

    with pytest.raises(units.DimensionalityError):
        z = mylist * units.ureg.hectare >= mylist * units.ureg.astronomical_unit


def test_addition_dimensionality_errors():
    with pytest.raises(units.DimensionalityError):
        x = 1.0 * units.angstrom + 1.0*units.ureg.kilograms

    with pytest.raises(units.DimensionalityError):
        y = 1.0 * units.fs - 1.0 * units.ureg.meter


def test_compatible_units_comparison():
    num = 1.0*units.angstrom
    assert abs(num-0.1*units.nm) < 1.0e-15 * units.angstrom
    assert 1.0 * units.ureg.meter > 123456.0 * units.nm


def test_default_units():
    assert units.default.length == units.angstrom
    assert units.default.mass == units.amu
    assert units.default.energy == units.eV


def test_default_unit_conversions():
    my_list = [1.0*units.angstrom, 1.0*units.nm, 1.0*units.a0]
    my_array = units.array(my_list).defunits()
    assert my_array.get_units() == units.default.convert(my_array).get_units()
    result = my_array.value_in(units.nm)
    np.testing.assert_almost_equal(result[0], 0.1, 9)
    np.testing.assert_almost_equal(result[1], 1.0, 9)
    np.testing.assert_almost_equal(result[2], 0.05291772, 7)


def test_scalar_comparisons_to_zero_ignore_units():
    num = 1.0*units.angstrom
    assert num > 0.0
    assert num > 0.0 * units.angstrom
    assert num > 0.0 * units.nm

    with pytest.raises(units.DimensionalityError):
        num > 0.0 * units.fs


def test_array_comparisons_to_zero_ignore_units():
    num = [1, -2, 0]*units.angstrom
    assert ((num > 0) == [True, False, False]).all()
    assert ((num == 0) == [False, False, True]).all()


def test_dimensionless_array_operations():
    arr = np.arange(5) * units.ureg.dimensionless
    assert (arr == [0, 1, 2, 3, 4]).all()

    arr[1:5] = [100, 101, 102, 103]
    assert (arr == [0, 100, 101, 102, 103]).all()

    assert arr[4].units == units.ureg.dimensionless
    assert arr[:3].units == units.ureg.dimensionless


def test_dimensionless_array_unit_checks():
    arr = np.arange(5) * units.ureg.dimensionless

    with pytest.raises(units.DimensionalityError):
        arr[0] = 5.0 * units.angstrom

    with pytest.raises(units.DimensionalityError):
        arr[:] = np.arange(5, 10) * units.angstrom

    arr[:] = np.arange(5, 10)
    assert (arr == np.arange(5, 10)).all()

    arr[:] = np.arange(10, 15) * units.ureg.dimensionless
    assert (arr == np.arange(10, 15)).all()


def test_array_unit_checks():
    arr = np.arange(5) * units.ureg.nm / units.ureg.fs

    with pytest.raises(units.DimensionalityError):
        arr[0] = 5.0 * units.angstrom

    with pytest.raises(units.DimensionalityError):
        arr[3] = 5.0

    with pytest.raises(units.DimensionalityError):
        arr[:] = np.arange(5, 10) * units.fs

    arr[2:3] = np.arange(5, 6) * units.ureg.angstrom / units.ureg.fs
    np.testing.assert_allclose(arr[2:3].magnitude,
                               np.arange(5, 6)/10.0)

    arr[:] = np.arange(10, 15) * units.ureg.micrometers / units.ureg.picoseconds
    np.testing.assert_allclose(arr.magnitude,
                               np.arange(10, 15))

def test_default_unit_conversions():
    assert abs(10.0 - (1.0*units.nm).defunits_value()) < 1e-10
    assert abs(1000.0 - (1.0*units.ps).defunits_value()) < 1e-10
    assert abs(1.0 - 6.022140857e23/((1.0*units.ureg.grams).defunits_value())) < 1e-6
    assert abs(103.642685491 - (1.0*units.angstrom**2*units.dalton/units.fs**2).defunits_value()
               ) < 1e-7
