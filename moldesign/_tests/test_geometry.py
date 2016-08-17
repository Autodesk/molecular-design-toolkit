""" Tests geometry routines
"""
import random

import itertools
import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u

registered_types = {}

# TODO: automated method testing based on its metadata - i.e. test to make sure parameters are
#       honored, test that it calcultes what it says it does, test that properties have the right
#       units and array shapes, etc.


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


def _make_mol_with_n_hydrogens(n):
    return mdt.Molecule([mdt.Atom('H') for i in xrange(n)])


def _apply_random_offsets(mol, idim):
    mol.positions[:, idim] += (random.random()-0.5)*100.0*u.angstrom


@typedfixture('atomcontainer')
def three_particle_right_angle():
    mol = _make_mol_with_n_hydrogens(3)
    mol.atoms[0].x = 1.0 * u.angstrom
    mol.atoms[2].y = 1.0 * u.angstrom

    for idim in xrange(3):
        _apply_random_offsets(mol, idim)

    return mol


@typedfixture('atomcontainer')
def four_particle_45_twist():
    mol = _make_mol_with_n_hydrogens(4)
    mol.positions = u.nm*[[0.1, 0.0, -0.5],
                          [0.0, 0.0, -0.5],
                          [0.0, 0.0, 0.5],
                          [0.2, -0.2, 0.5]]

    for idim in xrange(3):
        _apply_random_offsets(mol, idim)

    return mol


def test_dihedral_measure(four_particle_45_twist):
    mol = four_particle_45_twist
    np.testing.assert_almost_equal(mdt.dihedral(*mol.atoms).value_in(u.degrees),
                                   45.0,
                                   decimal=8)


def test_dihedral_sign_convention(four_particle_45_twist):
    mol = four_particle_45_twist
    mol.atoms[-1].y += 0.4 * u.nm
    np.testing.assert_almost_equal(mdt.dihedral(*mol.atoms).value_in(u.degrees),
                                   315.0,
                                   decimal=8)


def test_angle_measure(three_particle_right_angle):
    mol = three_particle_right_angle
    np.testing.assert_almost_equal(mdt.angle(*mol.atoms).value_in(u.degrees),
                                   90.0,
                                   decimal=8)


def test_distance_array(three_particle_right_angle):
    mol = three_particle_right_angle

    desired_distance_array = u.angstrom*[[0.0, 1.0, np.sqrt(2)],
                                         [1.0, 0.0, 1.0],
                                         [np.sqrt(2), 1.0, 0.0]]
    distance_array = mol.calc_distance_array()

    np.testing.assert_allclose(distance_array,
                               desired_distance_array,
                               atol=1e-8)


@pytest.mark.parametrize('objkey', registered_types['atomcontainer'])
def test_atomic_distance_measures_are_consistent(objkey, request):
    mol = request.getfuncargvalue(objkey)

    distance_array = mol.calc_distance_array()

    for i, j in itertools.product(xrange(3), xrange(3)):
        ai, aj = mol.atoms[i], mol.atoms[j]
        assert ai.distance(aj) == distance_array[i, j]
        assert mdt.distance(ai, aj) == distance_array[i, j]
        np.testing.assert_almost_equal(np.sum((ai.position - aj.position)**2).defunits_value(),
                                       (distance_array[i, j]**2).defunits_value(),
                                       decimal=10)
