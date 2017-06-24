import pytest

import moldesign as mdt
from moldesign import units as u
from .molecule_fixtures import h2
from .test_qm_xfaces import h2_rhfwfn

from .helpers import assert_almost_equal


@pytest.fixture
def h2_h2_with_ff(h2):
    ma = h2.copy()
    ma.residues[0].resname = 'UNL'
    ma.atoms[0].name = 'HA'
    ma.atoms[1].name = 'HB'
    list(ma.bonds)[0].align('x')

    mb = ma.copy()
    mb.translate([0.0, 2.0, 0.0]*u.angstrom)

    mol = ma.combine(mb)
    params = mdt.create_ff_parameters(ma, charges='gasteiger')
    params.assign(mol)
    return mol


def test_mechanical_embedding(h2_h2_with_ff):
    mol = h2_h2_with_ff

    mol.set_energy_model(mdt.models.MechanicalEmbeddingQMMM,
                         qm_atom_indices=[0, 1],
                         qm_model=mdt.models.RHF(basis='sto-3g'),
                         mm_model=mdt.models.OpenMMPotential)
    mol.calculate()
    qmprops = mol.properties.qmprops

    h2_qm = mdt.Molecule(mol.residues[0])
    h2_qm.set_energy_model(mdt.models.RHF, basis='sto-3g')
    h2_qm.calculate()

    assert abs(h2_qm.potential_energy - qmprops.potential_energy) < 1e-8 * u.hartree
    #assert_almost_equal(h2_qm.forces, qmprops.forces)
    assert_almost_equal(h2_qm.wfn.fock_ao, qmprops.wfn.fock_ao)
