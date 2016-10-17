import random

import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u

from . import helpers


registered_types = {}

def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@pytest.fixture
def small_molecule():
    return mdt.from_smiles('CNCOS(=O)C')


@typedfixture('mdready')
def parameterize_zeros(small_molecule):
    params = mdt.parameterize(small_molecule, charges='zero')
    mol = mdt.assign_forcefield(small_molecule, parameters=params)
    mol.set_energy_model(mdt.models.ForceField)
    return mol


@typedfixture('mdready')
def parameterize_am1bcc(small_molecule):
    params = mdt.parameterize(small_molecule, charges='am1-bcc', ffname='gaff')
    mol = mdt.assign_forcefield(small_molecule, parameters=params)
    mol.set_energy_model(mdt.models.ForceField)
    return mol


@typedfixture('mdready')
def protein_default_amber_forcefield():
    mol = mdt.from_pdb('1YU8')
    newmol = mdt.assign_forcefield(mol)
    newmol.set_energy_model(mdt.models.ForceField)
    return newmol


@typedfixture('mdready')
def gaff_model_gasteiger(small_molecule):
    small_molecule.set_energy_model(mdt.models.GAFF, charges='gasteiger')
    return small_molecule


@pytest.mark.parametrize('objkey', registered_types['mdready'])
def test_properties(objkey, request):
    mol = request.getfuncargvalue(objkey)
    energy = mol.calculate_potential_energy()
    forces = mol.calculate_forces()
    assert forces.shape == mol.positions.shape


@pytest.mark.skipif(mdt.interfaces.openmm.force_remote,
                    reason="Numerical derivatives need to be parallelized, "
                           "otherwise this takes too long")
@pytest.mark.parametrize('objkey', registered_types['mdready'])
def test_forces(objkey, request):
    mol = request.getfuncargvalue(objkey)

    anagrad = -mol.calculate_forces().defunits_value()
    numgrad = helpers.num_grad(mol,
                               mol.calculate_potential_energy,
                               step=0.005*u.angstrom
                               ).defunits_value()
    assert np.sqrt(np.sum((anagrad-numgrad) ** 2))/(3.0*mol.num_atoms) <= 1.0e-4  # this isn't good


@pytest.mark.parametrize('objkey', registered_types['mdready'])
def test_minimize(objkey, request):
    mol = request.getfuncargvalue(objkey)
    e1 = mol.calculate_potential_energy()
    mol = request.getfuncargvalue(objkey)
    traj = mol.minimize()
    assert mol.calculate_potential_energy() < e1
