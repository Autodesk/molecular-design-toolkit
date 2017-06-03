import random

import pytest

import moldesign as mdt
from moldesign import units as u

from . import helpers
from .molecule_fixtures import (pdb1yu8, small_molecule)
from .test_openbabel_xface import registered_types as obtypes
from .test_ambertools_xface import registered_types as ambtypes
from .test_ambertools_xface import (gaff_model_gasteiger, protein_default_amber_forcefield,
                                    parameterize_am1bcc, parameterize_zeros,
                                    protein_default_amber_forcefield)


registered_types = {}
registered_types.update(obtypes)
registered_types.update(ambtypes)


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@pytest.mark.parametrize('objkey', registered_types['hasmodel'])
def test_forces_and_energy_were_calculated(objkey, request):
    mol = request.getfixturevalue(objkey)
    energy = mol.calculate_potential_energy()
    forces = mol.calculate_forces()
    assert forces.shape == mol.positions.shape


@pytest.mark.skipif(mdt.interfaces.openmm.force_remote,
                    reason="Numerical derivatives need to be parallelized, "
                           "otherwise this takes too long")
@pytest.mark.parametrize('objkey', registered_types['hasmodel'])
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


@pytest.mark.parametrize('objkey', registered_types['hasmodel'])
def test_minimization_reduces_energy(objkey, request):
    mol = request.getfixturevalue(objkey)
    e1 = mol.calculate_potential_energy()
    traj = mol.minimize()
    assert mol.calculate_potential_energy() < e1
