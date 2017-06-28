import pytest
import moldesign as mdt

from . import helpers
from .molecule_fixtures import *


@pytest.mark.parametrize('fixturename', molecule_standards['hasmodel'])
def test_model_assigned(fixturename, request):
    mol = request.getfixturevalue(fixturename)
    assert mol.ff is not None


@pytest.mark.parametrize('objkey',
                         'ethylene ligand3aid mol_from_xyz mol_from_sdf'.split())
def test_parameterization_from_formats(objkey, request):
    mol = request.getfixturevalue(objkey)
    assert not mol.ff
    params = mdt.create_ff_parameters(mol, charges='gasteiger')
    assert params is not None
    _test_succesful_parameterization(mol)


@pytest.mark.screening
def test_parameterize_multiple_identical_small_molecules():
    m1 = mdt.from_smiles('O')
    params = mdt.create_ff_parameters(m1, charges='am1-bcc')
    assert params is not None
    m2 = m1.copy()
    m2.translate([4.0, 0.0, 0.0] * mdt.units.angstrom)
    mol = m1.combine(m2)
    ff = mdt.forcefields.DefaultAmber()
    ff.add_ff(params)

    ff.assign(mol)
    _test_succesful_parameterization(mol)


@pytest.mark.parametrize('chargemodel',
                         'esp gasteiger zero am1-bcc'.split())
def test_charge_models(ethylene, chargemodel):
    mol = ethylene
    if chargemodel == 'esp':
        pytest.xfail("ESP not yet implemented")
    assert not mol.ff
    params = mdt.create_ff_parameters(mol, charges=chargemodel)
    assert params is not None
    _test_succesful_parameterization(mol)


def _test_succesful_parameterization(mol):
    assert mol.ff
    mol.set_energy_model(mdt.models.ForceField)
    mol.calculate()
    assert 'potential_energy' in mol.properties
    assert 'forces' in mol.properties


def test_1yu8_default_amber_fix_and_assignment(protein_default_amber_forcefield):
    _test_succesful_parameterization(protein_default_amber_forcefield)


def test_ff_assignment_doesnt_change_topology(pdb3aid):
    m = pdb3aid
    protein = mdt.Molecule(m.get_atoms('protein'))
    ligand = mdt.Molecule(m.get_atoms('unknown'))
    ligff = mdt.create_ff_parameters(ligand, charges='gasteiger')

    mdt.guess_histidine_states(protein)
    mol = protein.combine(ligand)

    ff = mdt.forcefields.DefaultAmber()
    ff.add_ff(ligff)

    mdready = ff.create_prepped_molecule(mol)

    assert mdready.num_residues == mol.num_residues
    assert mdready.num_chains == mol.num_chains
    for c1, c2 in zip(mdready.chains, mol.chains):
        assert c1.name == c2.name
        assert c1.num_residues == c2.num_residues
        assert c1.index == c2.index

    for newr, oldr in zip(mdready.residues, mol.residues):
        assert newr.index == oldr.index
        if newr.resname == 'HIS':
            assert oldr.resname in 'HIS HID HIE HIP'.split()
        else:
            assert newr.resname == oldr.resname
        assert newr.pdbindex == oldr.pdbindex
        assert newr.chain.index == oldr.chain.index
        for atom in oldr:
            assert atom.name in newr
