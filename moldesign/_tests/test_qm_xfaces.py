""" Tests basic QM functionality and data structures
"""
import itertools
import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u

from . import helpers
from .molecule_fixtures import *

# TODO: automated method testing based on its metadata - i.e. test to make sure parameters are
#       honored, test that it calcultes what it says it does, test that properties have the right
#       units and array shapes, etc.


# Note that this is an incomplete set of models
models_to_test = list(itertools.product((mdt.models.NWChemQM, mdt.models.PySCFPotential),
                                        'sto-3g 6-31g'.split(),
                                        'rhf rks mp2'.split()))
model_ids = ['/'.join((model.__name__, theory, basis)) for (model, theory, basis) in models_to_test]

TESTSET = ['h2', 'small_molecule', 'benzene']


@pytest.fixture(params=models_to_test, ids=model_ids, scope='function')
def h2_with_model(request, h2):
    model, basis, theory = request.param

    if model is mdt.models.NWChemQM and theory == 'mp2':
        pytest.xfail('Not implemented')

    h2.set_energy_model(model, basis=basis, theory=theory)
    return h2


def test_minimization_trajectory(h2_with_model):
    mol = h2_with_model
    if mol.energy_model.params.theory == 'mp2':
        pytest.skip('Not testing mp2 minimizations at this time')

    assert 'potential_energy' not in mol.properties

    e1 = mol.calculate_potential_energy()
    p1 = mol.positions.copy()

    traj = mol.minimize()
    helpers.assert_something_resembling_minimization_happened(p1, e1, traj, mol)


@pytest.mark.parametrize('objkey', TESTSET)
def test_pyscf_rhf_sto3g_properties(objkey, request):
    mol = request.getfixturevalue(objkey)
    mol.set_energy_model(mdt.models.PySCFPotential, basis='sto-3g', theory='rhf')

    mol.calculate()

    assert 'potential_energy' in mol.properties
    assert 'wfn' in mol.properties
    assert 'canonical' in mol.wfn.orbitals
    assert 'atomic' in mol.wfn.orbitals
    assert mol.wfn.num_electrons == sum([atom.atnum for atom in mol.atoms]) \
                                    - mol.charge.value_in(u.q_e)


@pytest.mark.parametrize('objkey', TESTSET)
def test_pyscf_rhf_sto3g_matrices(objkey, request):
    mol = request.getfixturevalue(objkey)
    mol.set_energy_model(mdt.models.PySCFPotential, basis='sto-3g', theory='rhf')

    mol.calculate()
    basis = mol.wfn.aobasis
    canonical = mol.wfn.orbitals.canonical

    assert (mol.wfn.aobasis.fock == mol.wfn.fock_ao).all()
    assert (mol.wfn.orbitals.atomic.coeffs == np.identity(mol.wfn.nbasis)).all()

    np.testing.assert_allclose(canonical.to_ao(canonical.fock), mol.wfn.fock_ao, atol=1.e-9)
    np.testing.assert_allclose(canonical.from_ao(basis.overlaps), canonical.overlaps, atol=1.e-9)


def test_pyscf_casscf(h2):
    h2.set_energy_model(mdt.models.PySCFPotential,
                        theory='casscf',
                        active_electrons=2, active_orbitals=2,
                        state_average=3, basis='sto-3g')
    h2.calculate()
    # TODO: actually test results


def test_pyscf_casci(h2):
    pytest.xfail("CASCI interface doesn't yet work")
    h2.set_energy_model(mdt.models.PySCFPotential,
                        theory='casci',
                        active_electrons=2, active_orbitals=2, basis='6-31g')
    h2.calculate()
    # TODO: actually test results


@pytest.mark.parametrize('objkey', TESTSET)
def test_pyscf_rhf_sto3g_forces(objkey, request):
    mol = request.getfixturevalue(objkey)
    mol.set_energy_model(mdt.models.PySCFPotential, basis='sto-3g', theory='rhf')
    forces = mol.calc_forces()

    assert forces.shape == (mol.num_atoms, 3)


def test_calc_eri_tensor(h2):
    h2.set_energy_model(mdt.models.PySCFPotential, basis='sto-3g', theory='rhf')
    h2.calculate()
    eris = h2.wfn.aobasis.calc_eris()
    assert eris[0,0, 1,1] == eris[1,1, 0,0]
    assert eris[1,0, 0,1] == eris[0,1, 0,1]
    assert eris.nbasis == 2
    with pytest.raises(IndexError):
        eris[0,1,2,1]


@pytest.mark.screening
def test_aobasis(h2_rhf_sto3g):
    # it's sto-3g, so structure is simple
    aobasis = h2_rhf_sto3g.wfn.aobasis
    assert aobasis.basisname == 'sto-3g'
    assert (aobasis.coeffs == np.identity(2)).all()
    np.testing.assert_allclose(aobasis.fock, aobasis.fock.T)
    assert (aobasis.energies == aobasis.fock.diagonal()).all()

    assert aobasis.fock.dimensionality == u.eV.dimensionality

    for orb in h2_rhf_sto3g.wfn.aobasis:
        assert orb.aotype == '1s'
        assert orb.orbtype == 's'
        assert len(orb.primitives) == 3
        assert (orb.n, orb.m, orb.l) == (1, 0, 0)
        assert isinstance(orb.name, str)

# todo: deal with other shells, cartesian vs. spherical
