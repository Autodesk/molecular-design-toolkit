""" Tests the unit system
"""
import pytest
import numpy as np

from moldesign import units
from moldesign import units as u

# mark all tests in this module with these label (see ./conftest.py)
__PYTEST_MARK__ = ['internal', 'units']


@pytest.mark.screening
def test_scalar_comparison_dimensionality_errors():
    with pytest.raises(units.DimensionalityError):
        x = 1.0 * units.angstrom == 1.0*units.ureg.kilograms

    with pytest.raises(units.DimensionalityError):
        y = 1.0 * units.fs < 1.0 * units.ureg.meter

    with pytest.raises(units.DimensionalityError):
        z = 1.0 * units.ureg.hectare >= 1.0 * units.ureg.astronomical_unit


@pytest.mark.screening
def test_array_comparison_dimensionality_errors():
    mylist = [0.0, -1.0, 1.0]

    with pytest.raises(units.DimensionalityError):
        x = mylist * units.angstrom == mylist*units.ureg.kilograms

    with pytest.raises(units.DimensionalityError):
        y = mylist * units.fs < mylist * units.ureg.meter

    with pytest.raises(units.DimensionalityError):
        z = mylist * units.ureg.hectare >= mylist * units.ureg.astronomical_unit


@pytest.mark.screening
def test_addition_dimensionality_errors():
    with pytest.raises(units.DimensionalityError):
        x = 1.0 * units.angstrom + 1.0*units.ureg.kilograms

    with pytest.raises(units.DimensionalityError):
        y = 1.0 * units.fs - 1.0 * units.ureg.meter


@pytest.mark.screening
def test_compatible_units_comparison():
    num = 1.0*units.angstrom
    assert abs(num-0.1*units.nm) < 1.0e-15 * units.angstrom
    assert 1.0 * units.ureg.meter > 123456.0 * units.nm


def test_default_units():
    assert units.default.length == units.angstrom
    assert units.default.mass == units.amu
    assert units.default.energy == units.eV


@pytest.fixture
def si_units():
    system = units.UnitSystem(length=u.ureg.meters, mass=u.ureg.kg,
                              time=u.ureg.seconds, energy=u.ureg.joules)
    return system


@pytest.mark.parametrize(['dimension', 'expected', 'weird_units'],
                         [['length', u.ureg.meters, u.ureg.lightyear],
                          ['mass', u.ureg.kg, u.ureg.stone],
                          ['time', u.ureg.seconds, u.ureg.year],
                          ['momentum', u.ureg.kg * u.ureg.meter / u.ureg.second,
                                       u.ureg.lightyear * u.ureg.stone / u.ureg.year],
                          ['force', u.ureg.joule / u.ureg.meter, u.ureg.force_ton],
                          ['energy', u.ureg.joule, u.ureg.horsepower * u.ureg.seconds]],
                         ids='length mass time momentum force energy'.split())
def test_si_unit_system(si_units, dimension, expected, weird_units):
    assert getattr(si_units, dimension) == expected

    default_unit = getattr(u.default, dimension)
    my_quantity = 3.141 * default_unit
    si_quantity = si_units.convert(my_quantity)

    assert my_quantity.units != expected
    assert si_quantity.units == expected
    assert u.default.convert(abs(my_quantity - si_quantity)) < 1e-15 * default_unit

    setattr(si_units, dimension, weird_units)
    assert getattr(si_units, dimension) == weird_units

    weird_quantity = si_units.convert(my_quantity)
    assert abs(my_quantity-weird_quantity) < 1e-15 * default_unit
    assert weird_quantity.units == weird_units


@pytest.mark.screening
def test_default_unit_conversions():
    my_list = [1.0*units.angstrom, 1.0*units.nm, 1.0*units.a0]
    my_array = units.array(my_list).defunits()
    assert my_array.get_units() == units.default.convert(my_array).get_units()
    result = my_array.value_in(units.nm)
    np.testing.assert_almost_equal(result[0], 0.1, 9)
    np.testing.assert_almost_equal(result[1], 1.0, 9)
    np.testing.assert_almost_equal(result[2], 0.05291772, 7)


def test_dimensional_array_from_mixed_data():
    data = [np.arange(3) * u.angstrom,
            [1.0 * u.nm, 1.0 * u.nm**2/u.angstrom, 3.0*u.ureg.km],
            [1,1,1] * u.bohr]
    arr = u.array(data)

    bohr_in_ang = (1.0 * u.bohr).value_in(u.angstrom)
    expected = u.angstrom * [[0,1,2],
                             [10.0, 100.0, 3.0e13],
                             [bohr_in_ang, bohr_in_ang, bohr_in_ang]]
    np.testing.assert_allclose(arr, expected)
    assert isinstance(arr, units.MdtQuantity)


def test_dimensionless_array_from_mixed_data():
    data = [np.arange(3),
            [1.0 * u.dimensionless, 1.0 * u.nm/u.angstrom, 3.0],
            [1,1,1] * u.ureg.kg / u.ureg.g]
    arr = u.array(data)
    np.testing.assert_allclose(arr,
                               np.array([[0,1,2],
                                         [1,10,3],
                                         [1000,1000,1000]],dtype='float') )
    assert isinstance(arr, np.ndarray)


def test_inconsistent_array_units_raises_dimensionality_error():
    with pytest.raises(units.DimensionalityError):
        u.array([[1,2,3] * u.angstrom, [3*u.bohr, 4*u.eV, 5*u.bohr]])


def test_mixed_units_and_nonunits_raises_dimensionality_error():
    with pytest.raises(units.DimensionalityError):
        u.array([[1,2,3] * u.angstrom, [3*u.bohr, 4, 5*u.bohr]])


def test_no_nonsquare_arrays():
    with pytest.raises(ValueError):
        u.array([[1,2],[3]])


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


@pytest.mark.screening
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


@pytest.mark.parametrize(['a1','a2'],
                         ((np.ones(3), np.ones(3)*10/10),
                          (np.ones(3), u.array(np.ones(3))),
                          (np.ones(3)*u.nm, 10.0*np.ones(3)*u.angstrom),
                          (np.ones((3,3)) * u.nm / u.angstrom, np.ones(3) * 10)),
                         ids="numpy-numpy numpy-dimensionless nm-angstrom nm/angstrom-numpy".split()
                         )
def test_array_almost_equal_returns_true(a1, a2):
    assert u.arrays_almost_equal(a1, a2)
    assert u.arrays_almost_equal(a2, a1)


@pytest.mark.parametrize(['a1','a2'],
                         ((np.ones(3), np.ones(3)*10/10 * u.eV),
                          (np.ones(3)*u.nm, 10.0*np.ones(3)*u.dalton),
                          (u.array(np.ones(3)), np.ones(3) * u.kcalpermol)),
                         ids="numpy-ev nm-dalton dimensionless-kcalpermole".split()
                         )
def test_array_almost_equal_raises_dimensionality_error(a1, a2):
    with pytest.raises(u.DimensionalityError):
        u.arrays_almost_equal(a1, a2)
    with pytest.raises(u.DimensionalityError):
        u.arrays_almost_equal(a2, a1)



@pytest.mark.screening
def test_default_unit_conversions():
    assert abs(10.0 - (1.0*units.nm).defunits_value()) < 1e-10
    assert abs(1000.0 - (1.0*units.ps).defunits_value()) < 1e-10
    assert abs(1.0 - 6.022140857e23/((1.0*units.ureg.grams).defunits_value())) < 1e-6
    assert abs(103.642685491 - (1.0*units.angstrom**2*units.dalton/units.fs**2).defunits_value()
               ) < 1e-7


def test_getunits_doctests():
    assert units.get_units(1.0*units.angstrom) == units.MdtUnit('angstrom')
    assert units.get_units(np.array([1.0, 2, 3.0])) == units.MdtUnit('dimensionless')
    assert units.get_units([[1.0*units.dalton, 3.0*units.eV],
                            ['a'], 'gorilla']) == units.MdtUnit('amu')


@pytest.mark.screening
def test_setitem_from_quantity():
    myarray = np.arange(100) * u.angstrom
    with pytest.raises(u.DimensionalityError):
        myarray[:10] = np.arange(10)

    myarray[11] = 5.0 * u.ureg.meters

    myarray[:10] = np.arange(10) * u.nm
    for i in range(10):
        assert myarray[i].units == u.angstrom
        assert abs(myarray[i].magnitude - i*10) < 1e-14

    assert myarray[11].units == u.angstrom
    assert abs(myarray[11].magnitude - (5.0*u.ureg.meters).value_in(u.angstrom)) < 1e-14

    assert myarray.units == u.angstrom


@pytest.mark.screening
def test_setitem_from_list_of_quantities():
    myarray = np.arange(100) * u.angstrom
    myarray[:10] = list(np.arange(10) * u.nm)
    for i in range(10):
        assert myarray[i].units == u.angstrom
        assert abs(myarray[i].magnitude - i*10) < 1e-14 * u.angstrom


def test_setitem_slice_dimensionless():
    myarray = np.arange(100) * u.ureg.dimensionless
    assert isinstance(myarray, u.MdtQuantity)

    myarray[:10] = np.arange(10)
    assert myarray.units == u.ureg.dimensionless
    assert (myarray[:10] == np.arange(10)).all()


@pytest.fixture
def make_test_matrices():
    randommatrix = np.random.rand(10,7)
    randomvector = np.random.rand(10)
    unitmatrix = randommatrix * u.kcalpermol
    unitvector = randomvector * u.angstrom
    expected_matvec_product = randomvector.dot(randommatrix)
    return unitvector, unitmatrix, expected_matvec_product


def self_dot(a, b):
    return a.dot(b)


def self_ldot(a,b):
    return b.ldot(a)


def unit_dot(a,b):
    return u.dot(a, b)

MATRIX_MATHS = [self_dot, self_ldot, unit_dot]


@pytest.mark.parametrize('testfun', MATRIX_MATHS, ids=[x.__name__ for x in MATRIX_MATHS])
@pytest.mark.screening
def test_matrix_math_with_units(make_test_matrices, testfun):
    unitvector, unitmatrix, expected_matvec_product = make_test_matrices

    result = testfun(unitvector, unitmatrix)

    assert result.units == u.kcalpermol * u.angstrom
    np.testing.assert_array_almost_equal(expected_matvec_product,
                                         result.magnitude, decimal=12)

