import collections
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


@pytest.mark.parametrize('MinClass', (mdt.min.GradientDescent,
                                      mdt.min.BFGS,
                                      mdt.min.SmartMin))
def test_basic_minimization_remotely(harmonic_atom, MinClass):
    mol = harmonic_atom
    e0 = mol.calculate_potential_energy()
    p0 = mol.positions.copy()

    minimizer = MinClass(mol)
    traj = minimizer.runremotely()

    _check_basic_minimization_has_occured(p0, e0, traj)


MINIMIZERS = collections.OrderedDict([('gradient_descent', mdt.min.gradient_descent),
                                     ('leastsqr', mdt.min.sequential_least_squares),
                                     ('smart', mdt.min.minimize)])
# (ordered because pytest+xdist needs a definite ordering of parameters)


@pytest.mark.parametrize('minkey',(MINIMIZERS.keys()))
def test_constrained_distance_minimization(minkey):
    minimizer = MINIMIZERS[minkey]

    mol = mdt.Molecule([mdt.Atom(1), mdt.Atom(2)])
    mol.atoms[0].x = 2.0 * u.angstrom
    mol.atoms[1].x = 3.0 * u.angstrom
    mol.atoms[1].y = 2.0 * u.angstrom
    mol.set_energy_model(mdt.models.HarmonicOscillator, k=2.0*u.kcalpermol/u.angstrom**2)

    e0 = mol.calculate_potential_energy()
    p0 = mol.positions.copy()

    constraint = mol.constrain_distance(mol.atoms[0], mol.atoms[1])
    traj = minimizer(mol)

    _check_basic_minimization_has_occured(p0, e0, traj)
    assert abs(constraint.error()) < 1e-4 * u.angstrom


@pytest.mark.parametrize('minkey',(MINIMIZERS.keys()))
def test_constrained_dihedral_minimization(minkey):
    minimizer = MINIMIZERS[minkey]
    mol = mdt.from_smiles('C=C')
    assert mol.atoms[0].atnum == mol.atoms[1].atnum == 6  # make sure we're picking the right atoms

    dihedral = mdt.DihedralMonitor(mol.atoms[0], mol.atoms[1])
    assert dihedral.value == 0.0

    dihedral.value = 45 * u.degrees
    constraint = dihedral.constrain()

    mol.set_energy_model(mdt.models.OpenBabelPotential, forcefield='mmff94s')
    traj = minimizer(mol, nsteps=100)
    _check_basic_minimization_has_occured(traj.positions[0], traj.potential_energy[0], traj)

    assert constraint.error() <= 1.0 * u.degree
    traj_twists = traj.dihedral(mol.atoms[0], mol.atoms[1])
    assert (abs(traj_twists - 45 * u.degrees) <= 1.0 * u.degree).all()

    mol.clear_constraints()
    traj2 = minimizer(mol, nsteps=100)
    _check_basic_minimization_has_occured(traj2.positions[0], traj2.potential_energy[0], traj2)
    assert dihedral.value <= 5.0 * u.degrees


def _check_basic_minimization_has_occured(p0, e0, traj):
    assert traj.potential_energy[0] == e0
    assert traj.potential_energy[-1] < e0
    assert (traj.positions[0] == p0).all()
    assert (traj.positions[-1] != p0).any()
