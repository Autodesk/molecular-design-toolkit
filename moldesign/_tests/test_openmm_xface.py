import random
from itertools import product

import pytest

import moldesign as mdt
from moldesign import units as u

from . import helpers
from .molecule_fixtures import (pdb1yu8, small_molecule)
from .test_ambertools_xface import gaff_model_gasteiger, protein_default_amber_forcefield

# TODO: remove constraints from dynamics parameters - they should only live in the constraints array

TESTSYTEMS = ['small_mol', 'protein', 'protein_custom_constraints', 'protein_freeze_hbonds']
INTEGRATORS = ['verlet', 'langevin']


@pytest.fixture
def protein(protein_default_amber_forcefield):
    mol = protein_default_amber_forcefield
    mol.minimize(force_tolerance=0.5*u.eV/u.angstrom)  # perform a very partial minimization
    return mol


@pytest.fixture
def small_mol(gaff_model_gasteiger):
    mol = gaff_model_gasteiger
    mol.minimize(force_tolerance=0.5*u.eV/u.angstrom)  # perform a very partial minimization
    return mol


@pytest.fixture
def protein_custom_constraints(protein):
    # constrains distance between first and last c-alpha carbons plus position of c-alpha in PHE58
    mol = protein
    mol.constrain_distance(mol.chains['X'].n_terminal['CA'],
                           mol.chains['X'].c_terminal['CA'])
    mol.constrain_atom(mol.chains['X'].residues['PHE58']['CA'])
    return mol


@pytest.fixture
def protein_freeze_hbonds(protein):
    # constrains distance between first and last c-alpha carbons plus position of c-alpha in PHE58
    mol = protein
    mol.constrain_hbonds()
    return mol


@pytest.fixture
def langevin():
    return mdt.integrators.OpenMMLangevin(temperature=300.0*u.kelvin,
                                          constrain_hbonds=False,
                                          constrain_water=False)


@pytest.fixture
def verlet():
    return mdt.integrators.OpenMMVerlet(timestep=1.0 * u.fs,
                                        constrain_hbonds=False,
                                        constrain_water=False)


@pytest.mark.parametrize('objkey', TESTSYTEMS)
def test_forces_and_energy_were_calculated(objkey, request):
    mol = request.getfixturevalue(objkey)
    energy = mol.calculate_potential_energy()
    forces = mol.calculate_forces()
    assert forces.shape == mol.positions.shape


@pytest.mark.skipif(mdt.interfaces.openmm.force_remote,
                    reason="Numerical derivatives need to be parallelized, "
                           "otherwise this takes too long")
@pytest.mark.parametrize('objkey', TESTSYTEMS)
def test_analytical_vs_numerical_forces(objkey, request):
    mol = request.getfixturevalue(objkey)

    if mol.num_atoms > 20:
        atoms = random.sample(mol.atoms, 20)
    else:
        atoms = mol.atoms
    atom_indices = [atom.index for atom in atoms]

    anagrad = -mol.calculate_forces()[atom_indices]
    numgrad = helpers.num_grad(mol,
                               mol.calculate_potential_energy,
                               atoms=atoms,
                               step=0.005*u.angstrom)
    assert (anagrad-numgrad).norm()/(3.0*len(atoms)) <= 5.0e-4 * u.eV / u.angstrom


@pytest.mark.parametrize('objkey', TESTSYTEMS)
def test_minimization_reduces_energy(objkey, request):
    mol = request.getfixturevalue(objkey)
    p0 = mol.positions.copy()
    e0 = mol.calculate_potential_energy()
    traj = mol.minimize()
    if mol.constraints:
        pytest.xfail("This will fail until we replace minimizer with a CustomIntegrator")
    helpers.assert_something_resembling_minimization_happened(p0, e0, traj, mol)


@pytest.mark.parametrize('systemkey,integratorkey', product(TESTSYTEMS, INTEGRATORS))
def test_openmm_dynamics(systemkey, integratorkey, request):
    mol = request.getfixturevalue(systemkey)
    mol.set_integrator(request.getfixturevalue(integratorkey))

    mol.integrator.prep()
    assert mol.integrator.sim.system is mol.energy_model.sim.system
    assert mol.integrator.sim.system.getNumConstraints() == len(
            [c for c in mol.constraints if c.desc == 'position'])

    p0 = mol.positions.copy()
    t0 = mol.time
    traj = mol.run(10.0 * u.ps)

    # Basic dynamics and constraint sanity checks:
    helpers.assert_something_resembling_dynamics_happened(traj, mol, p0, t0, 10.0*u.ps)

    if 'temperature' in mol.integrator.params:
        temp = mol.integrator.params.temperature
        if mol.num_atoms > 50:
            assert (temp / 2.0 <= traj.kinetic_temperature[-1] <= 500 * u.kelvin)
        else:  # small molecules have a big range, just a very loose sanity check here
            assert mol.kinetic_temperature > temp / 10.0

    else:  # Assume constant energy dynamics (may need to check explicitly in future)
        energy = traj.potential_energy+traj.kinetic_energy
        assert abs(energy[0]-energy[-1]) <= 0.01*mol.num_atoms*u.kcalpermol


def test_cleared_constraints_are_no_longer_applied(protein_custom_constraints, langevin):
    mol = protein_custom_constraints

    t0 = mol.time
    p0 = mol.positions.copy()
    mol.set_integrator(langevin)
    traj = mol.run(2.0 * u.ps)

    helpers.assert_something_resembling_dynamics_happened(traj, mol, p0, t0, 2*u.ps)

    for constraint in mol.constraints:
        assert constraint.satisfied()

    oldconstraints = mol.constraints[:]

    assert len(mol.constraints) > 0
    mol.clear_constraints()
    assert len(mol.constraints) == 0

    t1 = mol.time
    p1 = mol.positions.copy()
    traj = mol.run(10.0 * u.ps)

    for constraint in oldconstraints:
        assert not constraint.satisfied()  # it would be very very unlikely, I think

    helpers.assert_something_resembling_dynamics_happened(traj, mol, p1, t1, 5*u.ps)


@pytest.mark.parametrize('integkey', INTEGRATORS)
def test_unsupported_constraint_types(protein, integkey, request):
    integrator = request.getfixturevalue(integkey)
    protein.set_integrator(integrator)

    assert protein.energy_model._prepped
    protein.constrain_dihedral(*protein.atoms[:4])
    assert not protein.energy_model._prepped

    protein.calculate(use_cache=False)  # this should work, because there's no motion invovled

    with pytest.raises(mdt.exceptions.NotSupportedError):
        traj = protein.run(1*u.ps)

    t0 = protein.time
    p0 = protein.positions.copy()
    protein.clear_constraints()
    traj = protein.run(1*u.ps)  # should NOT raise a fuss now
    helpers.assert_something_resembling_dynamics_happened(traj, protein, p0, t0, 1*u.ps)

    protein.constrain_angle(*protein.atoms[:3])
    with pytest.raises(mdt.exceptions.NotSupportedError):
        protein.run(1*u.ps)


@pytest.mark.parametrize('integkey', INTEGRATORS)
def test_fallback_to_builtin_minimizer_for_arbitrary_constraints(small_mol, integkey, request):
    mol = small_mol
    assert len(mol.constraints) == 0
    mol.set_integrator(request.getfixturevalue(integkey))

    mol.constrain_angle(*mol.atoms[:3])
    assert len(mol.constraints) == 1

    p0 = mol.positions.copy()
    e0 = mol.calculate_potential_energy()

    traj = mol.minimize()
    helpers.assert_something_resembling_minimization_happened(p0, e0, traj, mol)


def test_list_platforms():  # doesn't do much right now
    platforms = mdt.interfaces.openmm.list_openmmplatforms()
    print('Found platforms %d: ', (len(platforms), platforms))
