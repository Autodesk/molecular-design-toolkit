import pytest

import moldesign as mdt
from moldesign import units as u

from . import helpers

# Tests:
#  1. internal bonds on QM region are removed in all cases
#  2. wavefunction is perturbed for electrostatic embedding
#  3.



@pytest.fixture
def h2params():
    mol = mdt.from_smiles('[H][H]')
    mol.atoms[0].name = 'HA'
    mol.atoms[1].name = 'HB'
    mol.residues[0].name = 'H2'

    params = mdt.create_ff_parameters(mol, charges='gasteiger')
    return mol, params


@pytest.fixture(scope='function')
def h2_h2_with_ff(h2params):
    ma, params = h2params
    ma.residues[0].resname = 'UNL'
    ma.atoms[0].name = 'HA'
    ma.atoms[1].name = 'HB'
    list(ma.bonds)[0].align('x')

    mb = ma.copy()
    mb.translate([0.0, 2.0, 0.0]*u.angstrom)

    mol = ma.combine(mb)
    params.assign(mol)
    return mol


@pytest.fixture
def h2_mm(h2params):
    mol, params = h2params
    params.assign(mol)
    mol.set_energy_model(mdt.models.OpenMMPotential)
    return mol


@pytest.fixture
def h2_qm(h2params):
    mol, params = h2params
    mol.set_energy_model(mdt.models.RHF, basis='sto-3g')
    return mol


@pytest.fixture
def h2_h2_mm(h2_h2_with_ff):
    h2_h2_with_ff.set_energy_model(mdt.models.OpenMMPotential)
    return h2_h2_with_ff


@pytest.fixture
def h2_h2_rhf(h2_h2_with_ff):
    h2_h2_with_ff.set_energy_model(mdt.models.RHF, basis='sto-3g')
    return h2_h2_with_ff


@pytest.fixture
def h2_h2_mechanical_embedding_rhf(h2_h2_with_ff):
    mol = h2_h2_with_ff
    mol.set_energy_model(mdt.models.MechanicalEmbeddingQMMM,
                         qm_atom_indices=[0, 1],
                         qm_model=mdt.models.RHF(basis='sto-3g'),
                         mm_model=mdt.models.OpenMMPotential)
    return mol


@pytest.fixture
def h2_h2_mechanical_embedding_zeroqm(h2_h2_with_ff):
    mol = h2_h2_with_ff
    mol.set_energy_model(mdt.models.MechanicalEmbeddingQMMM,
                         qm_atom_indices=[0, 1],
                         qm_model=helpers.ZeroEnergy,
                         mm_model=mdt.models.OpenMMPotential)
    return mol


def test_mechanical_embedding_wfn(h2_h2_mechanical_embedding_rhf):
    mol = h2_h2_mechanical_embedding_rhf

    mol.calculate()
    qmprops = mol.properties.qmprops
    mmprops = mol.properties.mmprops

    h2_qm = mdt.Molecule(mol.residues[0])
    h2_qm.set_energy_model(mdt.models.RHF, basis='sto-3g')
    h2_qm.calculate()

    assert abs(h2_qm.potential_energy - qmprops.potential_energy) < 1e-8 * u.hartree
    helpers.assert_almost_equal(h2_qm.wfn.fock_ao, qmprops.wfn.fock_ao)
    assert qmprops.potential_energy + mmprops.potential_energy == mol.potential_energy
