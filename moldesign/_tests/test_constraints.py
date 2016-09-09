""" Tests constraint routines
"""
import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u


def _make_mol_with_n_hydrogens(n):
    return mdt.Molecule([mdt.Atom('H') for i in xrange(n)])

@pytest.fixture
def four_particle_45_twist():
    mol = _make_mol_with_n_hydrogens(4)
    mol.positions = u.nm*[[0.1, 0.0, -0.5],
                          [0.0, 0.0, -0.5],
                          [0.0, 0.0, 0.5],
                          [0.2, -0.2, 0.5]]

    for iatom in xrange(3):
        mol.atoms[iatom].bond_to(mol.atoms[iatom+1], 1)

    return mol


def test_dihedral_constraint_errors(four_particle_45_twist):
    mol = four_particle_45_twist
    constraint = mol.constrain_dihedral(*mol.atoms)

    assert constraint.error() == 0.0

    constraint.value = 30.0 * u.degrees
    np.testing.assert_allclose(constraint.error().value_in(u.degrees), 15.0)

    constraint.value = 60.0 * u.degrees
    np.testing.assert_allclose(constraint.error().value_in(u.degrees), -15.0)


def test_dihedral_constraint_errors_at_0(four_particle_45_twist):
    mol = four_particle_45_twist
    mol.atoms[3].position = [0.1, 0.0, 0.5] * u.angstrom

    constraint = mol.constrain_dihedral(*mol.atoms)
    assert constraint.error() == 0.0

    for angle in (0, -10, 10, 90, -90) * u.degrees:
        constraint.value = angle
        np.testing.assert_allclose(constraint.error().value_in(u.degrees), -angle.magnitude)
