""" Tests geometry routines
"""
from builtins import range
import random

import itertools
import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u

from . import helpers

registered_types = {}


__PYTEST_MARK__ = 'internal'  # mark all tests in this module with this label (see ./conftest.py)


# TODO: automated method testing based on its metadata - i.e. test to make sure parameters are
#       honored, test that it calcultes what it says it does, test that properties have the right
#       units and array shapes, etc.

# step for numerical derivative testing


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


def _make_mol_with_n_hydrogens(n):
    return mdt.Molecule([mdt.Atom('H') for i in range(n)])


def _apply_random_offsets(mol, idim):
    mol.positions[:, idim] += (random.random()-0.5)*100.0*u.angstrom



@typedfixture('atomcontainer', scope='function')
def three_particle_right_angle():
    mol = _make_mol_with_n_hydrogens(3)
    mol.atoms[0].x = 1.0 * u.angstrom
    mol.atoms[2].y = 1.0 * u.angstrom

    for idim in range(3):
        _apply_random_offsets(mol, idim)

    return mol


@typedfixture('atomcontainer', scope='function')
def four_particle_45_twist():
    mol = _make_mol_with_n_hydrogens(4)
    mol.positions = u.nm*[[0.1, 0.0, -0.5],
                          [0.0, 0.0, -0.5],
                          [0.0, 0.0, 0.5],
                          [0.2, -0.2, 0.5]]

    for idim in range(3):
        _apply_random_offsets(mol, idim)

    for iatom in range(3):
        mol.atoms[iatom].bond_to(mol.atoms[iatom+1], 1)

    return mol

########################
# Dihedrals            #
########################
def test_dihedral_measure(four_particle_45_twist):
    mol = four_particle_45_twist
    np.testing.assert_almost_equal(mdt.dihedral(*mol.atoms).value_in(u.degrees),
                                   45.0,
                                   decimal=8)


def test_dihedral_two_atom_selection(four_particle_45_twist):
    mol = four_particle_45_twist
    np.testing.assert_almost_equal(mdt.dihedral(*mol.atoms[1:3]).value_in(u.degrees),
                                   45.0,
                                   decimal=8)

    with pytest.raises(ValueError):  # raises exception because it's not part of a dihedral
        mdt.dihedral(mol.atoms[0], mol.atoms[1])


def test_set_dihedral(four_particle_45_twist):
    mol = four_particle_45_twist
    mdt.set_dihedral(mol.atoms[0], mol.atoms[1], mol.atoms[2], mol.atoms[3], 10.0 * u.degrees)
    np.testing.assert_almost_equal(mdt.dihedral(*mol.atoms).value_in(u.degrees),
                                   10.0,
                                   decimal=8)


@pytest.mark.screening
def test_set_dihedral_sign_convention(four_particle_45_twist):
    mol = four_particle_45_twist
    mdt.set_dihedral(mol.atoms[0], mol.atoms[1], mol.atoms[2], mol.atoms[3], -23.0 * u.degrees)
    np.testing.assert_almost_equal(mdt.dihedral(*mol.atoms).value_in(u.degrees),
                                   337.0,
                                   decimal=8)


def test_set_dihedral_two_atom_selection(four_particle_45_twist):
    mol = four_particle_45_twist
    mdt.set_dihedral(mol.atoms[1], mol.atoms[2], 10.0 * u.degrees)
    np.testing.assert_almost_equal(mdt.dihedral(*mol.atoms).value_in(u.degrees),
                                   10.0,
                                   decimal=8)

    with pytest.raises(ValueError):  # raises exception because it's not part of a dihedral
        mdt.set_dihedral(mol.atoms[0], mol.atoms[1], 5.0 * u.degrees)


def test_set_dihedral_bond_no_adjust(four_particle_45_twist):
    mol = four_particle_45_twist
    bond = mdt.Bond(mol.atoms[1], mol.atoms[2])
    mdt.set_dihedral(bond, 10.0 * u.degrees, adjustmol=False)
    np.testing.assert_almost_equal(mdt.dihedral(*mol.atoms).value_in(u.degrees),
                                   10.0,
                                   decimal=8)


def test_dihedral_sign_convention(four_particle_45_twist):
    mol = four_particle_45_twist
    mol.atoms[-1].y += 0.4 * u.nm
    np.testing.assert_almost_equal(mdt.dihedral(*mol.atoms).value_in(u.degrees),
                                   315.0,
                                   decimal=8)

# TODO: test behavior at discontinuities (180, -180)

@pytest.mark.screening
def test_dihedral_gradient(four_particle_45_twist):
    mol = four_particle_45_twist

    dihe = mdt.DihedralMonitor(*mol.atoms)
    calc_grad = dihe.gradient()
    num_grad = helpers.num_grad(mol, lambda: dihe.value)

    np.testing.assert_allclose(calc_grad.defunits_value(),
                               num_grad.defunits_value(),
                               atol=5.0*helpers.DEFSTEP.defunits_value())


def test_dihedral_gradient_sign_convention(four_particle_45_twist):
    mol = four_particle_45_twist
    mol.atoms[-1].y += 0.4 * u.nm
    dihe = mdt.DihedralMonitor(*mol.atoms)
    calc_grad = dihe.gradient()
    num_grad = helpers.num_grad(mol, lambda: dihe.value)

    np.testing.assert_allclose(calc_grad,
                               num_grad,
                               atol=5.0*helpers.DEFSTEP.defunits_value())


########################
# Angles               #
########################
def test_angle_measure(three_particle_right_angle):
    mol = three_particle_right_angle
    np.testing.assert_almost_equal(mdt.angle(*mol.atoms).value_in(u.degrees),
                                   90.0,
                                   decimal=8)


@pytest.mark.screening
def test_angle_gradient(three_particle_right_angle):
    mol = three_particle_right_angle

    ang = mdt.AngleMonitor(*mol.atoms)
    assert abs(ang.value.value_in(u.degrees) - 90.0) <= 1.0e-8
    calc_grad = ang.gradient()
    num_grad = helpers.num_grad(mol, lambda:ang.value)

    np.testing.assert_allclose(calc_grad.defunits_value(),
                               num_grad.defunits_value(),
                               atol=5.0*helpers.DEFSTEP.defunits_value())


def test_set_angle_with_monitor(three_particle_right_angle):
    mol = three_particle_right_angle
    ang = mdt.AngleMonitor(*mol.atoms)

    ang.value = 45 * u.degrees
    assert abs(mdt.angle(*mol.atoms) - (45 * u.degrees)) < 0.1 * u.degrees


def test_set_angle_noadjust(four_particle_45_twist):
    mol = four_particle_45_twist
    assert mdt.angle(*mol.atoms[:3]) == 90.0 * u.degrees
    final = 45 * u.degrees
    origpos = mol.positions.copy()
    mdt.set_angle(mol.atoms[0], mol.atoms[1], mol.atoms[2], final, adjustmol=False)

    assert abs(mdt.angle(*mol.atoms[:3]) - final) < 0.1 * u.degrees
    assert (mol.positions[-1] == origpos[-1]).all()



########################
# Distances            #
########################
def test_distance_array(three_particle_right_angle):
    mol = three_particle_right_angle

    desired_distance_array = u.angstrom*[[0.0, 1.0, np.sqrt(2)],
                                         [1.0, 0.0, 1.0],
                                         [np.sqrt(2), 1.0, 0.0]]
    distance_array = mol.calc_distance_array()

    np.testing.assert_allclose(distance_array,
                               desired_distance_array,
                               atol=1e-8)


def test_set_distance_and_adjust(four_particle_45_twist):
    mol = four_particle_45_twist
    origpos = mol.positions.copy()
    distance = mdt.DistanceMonitor(mol.atoms[1], mol.atoms[2])
    olddist = distance.value
    distance.value *= 2.0

    displacement = np.sqrt(((origpos[0] - mol.positions[0])**2).sum()) + \
                   np.sqrt(((origpos[-1] - mol.positions[-1])**2).sum())
    assert abs(mdt.distance(mol.atoms[1], mol.atoms[2]) - 2.0*olddist) <= 1e-9 * u.angstrom
    assert abs(displacement - olddist) < 1.0e-9 * u.angstrom


def test_set_distance_noadjust(four_particle_45_twist):
    mol = four_particle_45_twist
    origpos = mol.positions.copy()
    olddist = mdt.distance(mol.atoms[1], mol.atoms[2])
    mdt.set_distance(mol.atoms[1], mol.atoms[2], 2.0 * olddist, adjustmol=False)

    assert abs(mdt.distance(mol.atoms[1], mol.atoms[2]) - 2.0*olddist) <= 1e-9 * u.angstrom
    assert (origpos[0] == mol.positions[0]).all() and (origpos[-1] == mol.positions[-1]).all()


@pytest.mark.parametrize('objkey', registered_types['atomcontainer'])
def test_atomic_distance_measures_are_consistent(objkey, request):
    mol = request.getfixturevalue(objkey)

    distance_array = mol.calc_distance_array()

    for i, j in itertools.product(range(3), range(3)):
        ai, aj = mol.atoms[i], mol.atoms[j]
        assert ai.distance(aj) == distance_array[i, j]
        assert mdt.distance(ai, aj) == distance_array[i, j]
        np.testing.assert_almost_equal(np.sum((ai.position - aj.position)**2).defunits_value(),
                                       (distance_array[i, j]**2).defunits_value(),
                                       decimal=10)


def test_center_of_mass(four_particle_45_twist):
    mol = four_particle_45_twist

    mol.positions = u.nm*[[0.1, 0.0, -0.5],
                          [0.0, 0.0, -0.5],
                          [0.0, 0.0, 0.5],
                          [0.2, -0.2, 0.5]]

    desired_com_angstroms = np.array([0.1+0.2, -0.2, 0.0]) * 10.0 / 4.0
    np.testing.assert_almost_equal(mol.center_of_mass.defunits_value(),
                                   desired_com_angstroms)

    mol.atoms[0].mass = 5.0 * u.ureg.kilograms
    mol.atoms[1].mass = 10.0 * u.ureg.kilograms
    mol.atoms[2].mass = 5.0 * u.ureg.kilograms
    mol.atoms[3].mass = 10.0 * u.ureg.kilograms

    desired_com_angstroms = np.array([0.1+0.4, -0.4, 0.0]) * 10.0 / 6.0
    np.testing.assert_almost_equal(mol.center_of_mass.defunits_value(),
                                   desired_com_angstroms)


def test_set_center_of_mass(four_particle_45_twist):
    # reset COM
    four_particle_45_twist.com = [0, 0, 0] * u.angstrom
    np.testing.assert_almost_equal(four_particle_45_twist.com.value_in(u.angstrom),
                                   ([0, 0, 0] * u.angstrom).value_in(u.angstrom))

    # set it to its current position
    four_particle_45_twist.com = [0, 0, 0] * u.angstrom
    np.testing.assert_almost_equal(four_particle_45_twist.com.value_in(u.angstrom),
                                   ([0, 0, 0] * u.angstrom).value_in(u.angstrom))

    # move COM elsewhere
    four_particle_45_twist.com = [10, 0, -10] * u.angstrom
    np.testing.assert_almost_equal(four_particle_45_twist.com.value_in(u.angstrom),
                                   ([10, 0, -10] * u.angstrom).value_in(u.angstrom))


def test_distance_gradient(three_particle_right_angle):
    mol = three_particle_right_angle

    dist = mdt.DistanceMonitor(*mol.atoms[:2])
    assert dist.value == mol.atoms[0].distance(mol.atoms[1])

    calc_grad = dist.gradient()
    num_grad = helpers.num_grad(mol, lambda:dist.value)

    np.testing.assert_allclose(calc_grad.defunits_value(),
                               num_grad.defunits_value(),
                               atol=5.0*helpers.DEFSTEP.defunits_value())

#########################
# Utilities             #
#########################
def test_sub_angles():
    from moldesign.geom import sub_angles
    np.testing.assert_allclose(sub_angles(np.pi*u.radian, np.pi/2.0*u.radian),
                               np.pi/2.0 * u.radian)

    np.testing.assert_allclose(sub_angles(360*u.degrees, 179*u.degrees),
                               -179*u.degrees)

    np.testing.assert_allclose(sub_angles(360*u.degrees, 3*u.degrees),
                               -3*u.degrees)

    np.testing.assert_allclose(sub_angles(720*u.degrees, -360*u.degrees),
                               0*u.degrees)

    np.testing.assert_allclose(sub_angles(180*u.degrees, 270*u.degrees),
                               -90*u.degrees)

    np.testing.assert_allclose(sub_angles(270*u.degrees, 0*u.degrees),
                               -90*u.degrees)
