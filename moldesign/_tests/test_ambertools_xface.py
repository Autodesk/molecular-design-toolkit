import pytest
import moldesign as mdt

from . import helpers
from .molecule_fixtures import *

registered_types = {}


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


# TODO: tests for parameterization error detection


@pytest.fixture(scope='function')
def mol_from_pdb():
    mol = mdt.from_pdb('3aid')
    return mdt.Molecule(mol.chains['A'].residues['ARQ401'])


@pytest.fixture(scope='function')
def mol_from_smiles():
    return mdt.from_smiles('CC')


@pytest.mark.parametrize('objkey',
                         'ethylene ligand3aid mol_from_xyz mol_from_sdf'.split())
def test_parameterization_from_formats(objkey, request):
    mol = request.getfixturevalue(objkey)
    assert not mol.ff
    params = mdt.create_ff_parameters(mol, charges='gasteiger')
    assert params is not None
    _test_succesful_parameterization(mol)


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
def test_charge_models(mol_from_smiles, chargemodel):
    mol = mol_from_smiles
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


def test_1yu8_default_amber_fix_and_assignment(pdb1yu8):
    mol = pdb1yu8
    ff = mdt.forcefields.DefaultAmber()
    newmol = ff.create_prepped_molecule(mol)
    _test_succesful_parameterization(newmol)


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



@typedfixture('hasmodel')
def parameterize_zeros(small_molecule):
    return _param_small_mol(small_molecule, 'zero')


@typedfixture('hasmodel')
def parameterize_am1bcc(small_molecule):
    return _param_small_mol(small_molecule, 'am1-bcc')


@typedfixture('hasmodel')
def gaff_model_gasteiger(small_molecule):
    small_molecule.set_energy_model(mdt.models.GaffSmallMolecule, partial_charges='gasteiger')
    small_molecule.energy_model.prep()
    return small_molecule


def _param_small_mol(small_molecule, chargemodel):
    params = mdt.create_ff_parameters(small_molecule, charges=chargemodel, baseff='gaff2')
    params.assign(small_molecule)
    small_molecule.set_energy_model(mdt.models.ForceField)
    return small_molecule


@typedfixture('hasmodel')
def protein_default_amber_forcefield(pdb1yu8):
    mol = pdb1yu8
    ff = mdt.forcefields.DefaultAmber()
    newmol = ff.create_prepped_molecule(mol)
    newmol.set_energy_model(mdt.models.ForceField)
    return newmol


@pytest.mark.parametrize('fixturename', registered_types['hasmodel'])
def test_model_assigned(fixturename, request):
    mol = request.getfixturevalue(fixturename)
    assert mol.ff is not None
