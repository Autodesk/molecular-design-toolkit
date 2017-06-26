import random
from itertools import product

import pytest

import moldesign as mdt
from moldesign import units as u

from . import helpers
from .molecule_fixtures import (pdb1yu8, small_molecule)
from .test_ambertools_xface import gaff_model_gasteiger, protein_default_amber_forcefield

# TODO: remove constraints from dynamics parameters - they should only live in the constraints array

TESTSYTEMS = ['gaff_model_gasteiger', 'protein_default_amber_forcefield',
              'protein_custom_constraints', 'protein_freeze_hbonds']
INTEGRATORS = ['verlet', 'langevin']


@pytest.fixture
def protein_custom_constraints(protein_default_amber_forcefield):
    # constrains distance between first and last c-alpha carbons plus position of c-alpha in PHE58
    mol = protein_default_amber_forcefield
    mol.constrain_distance(mol.chains['X'].n_terminal['CA'],
                           mol.chains['X'].c_terminal['CA'])
    mol.constrain_atom(mol.chains['X'].residues['PHE58']['CA'])
    return mol


@pytest.fixture
def protein_freeze_hbonds(protein_default_amber_forcefield):
    # constrains distance between first and last c-alpha carbons plus position of c-alpha in PHE58
    mol = protein_default_amber_forcefield
    mol.constrain_hbonds()
    return mol


@pytest.fixture
def langevin():
    return mdt.integrators.OpenMMLangevin(temperature=300.0*u.kelvin, constrain_hbonds=False)


@pytest.fixture
def verlet():
    return mdt.integrators.OpenMMVerlet(timestep=1.0 * u.fs, constrain_hbonds=False)


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
    mol.set_integrator(langevin)
    traj = mol.run(10.0 * u.ps)

    for constraint in mol.constraints:
        assert constraint.satisfied()

    constraint_refs = mol.constraints[:]

    mol.clear_constraints()
    traj = mol.run(20.0 * u.ps)

    for constraint in mol.constraints:
        assert not constraint.satisfied()  # it would be very very unlikely, I think


def test_unsupported_constraint_types(protein):
    protein.constrain_dihedral(*protein.atoms[:4])
    with pytest.raises(mdt.exceptions.NotSupportedError):
        protein.calculate()

    protein.clear_constraints()
    protein.calculate()  # should NOT raise a fuss now

    protein.constrain_angle(*protein.atoms[:3])
    with pytest.raises(mdt.exceptions.NotSupportedError):
        protein.calculate()


def test_list_platforms():  # doesn't do much right now
    platforms = mdt.interfaces.openmm.list_openmmplatforms()
    print('Found platforms %d: ', (len(platforms), platforms))
