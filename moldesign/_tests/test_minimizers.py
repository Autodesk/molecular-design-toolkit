import pytest


import moldesign as mdt
from moldesign import units as u


@pytest.fixture(scope='function')
def harmonic_atom():
    mol = mdt.Molecule([mdt.Atom(1)])
    mol.atoms[0].x = 2.0 * u.angstrom
    mol.set_energy_model(mdt.models.HarmonicOscillator, k=2.0*u.kcalpermol/u.angstrom**2)
    return mol


@pytest.mark.parametrize('MinClass', (mdt.min.GradientDescent,
                                      mdt.min.BFGS,
                                      mdt.min.SmartMin))
def test_basic_minimization(harmonic_atom, MinClass):
    mol = harmonic_atom
    e0 = mol.calculate_potential_energy()
    p0 = mol.positions.copy()

    minimizer = MinClass(mol)
    traj = minimizer()

    _check_basic_minimization_has_occured(p0, e0, traj)


def test_constrained_minimization():
    mol = mdt.Molecule([mdt.Atom(1), mdt.Atom(2)])
    mol.atoms[0].x = 2.0 * u.angstrom
    mol.atoms[1].x = 3.0 * u.angstrom
    mol.atoms[1].y = 2.0 * u.angstrom
    mol.set_energy_model(mdt.models.HarmonicOscillator, k=2.0*u.kcalpermol/u.angstrom**2)

    e0 = mol.calculate_potential_energy()
    p0 = mol.positions.copy()

    constraint = mol.constrain_distance(mol.atoms[0], mol.atoms[1])
    traj = mol.minimize()  # use the built-in logic for this one

    _check_basic_minimization_has_occured(p0, e0, traj)
    assert abs(constraint.error()) < 1e-4 * u.angstrom


def _check_basic_minimization_has_occured(p0, e0, traj):
    assert traj.potential_energy[0] == e0
    assert traj.potential_energy[-1] < e0
    assert (traj.positions[0] == p0).all()
    assert (traj.positions[-1] != p0).any()
