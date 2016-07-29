import pytest

from moldesign import units

u = units  # because I'm lazy -- AMV


def test_comparisons_to_zero():
    num = 1.0*units.angstrom
    assert num > 0.0
    assert num > 0.0 * units.angstrom
    assert num > 0.0 * units.nm
    with pytest.raises(units.DimensionalityError):
        tmp = 1.0*u.angstrom == 1.0 * u.ureg.kilograms


def test_compatible_units_comparison():
    num = 1.0*units.angstrom
    assert abs(num-0.1*units.nm) < 1.0e-15 * u.angstrom
    assert 1.0 * units.ureg.meter > 123456.0 * units.nm

