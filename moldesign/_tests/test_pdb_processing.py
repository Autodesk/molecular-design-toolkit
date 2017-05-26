"""
This file collects tests relating to the many special cases we deal with when
trying to process biomolecular structures from the pdb.

These currently just test the specific cases that we've implemented fixes for
"""
import pytest
import moldesign as mdt
from .helpers import get_data_path

import parmed

from distutils.version import LooseVersion


@pytest.fixture
def pdb_3ac2():
    return mdt.read(get_data_path('3ac2.pdb'))


@pytest.fixture
def pdb_3ac2_roundtrip(pdb_3ac2):
    return mdt.read(pdb_3ac2.write('pdb'), format='pdb')


@pytest.mark.parametrize('mol', 'pdb_3ac2 pdb_3ac2_roundtrip'.split())
def test_missing_terminal_atoms_3ac2(request, mol):
    """ Tests that we can still detect terminal residues even if the peptide-bonded atoms
    are missing from the structure
    """
    mol = request.getfixturevalue(mol)
    assert mol.chains['A'].n_terminal is not None
    assert mol.chains['A'].c_terminal is not None


@pytest.fixture
def pdb_1hpk():
    return mdt.read(get_data_path('1hpk.pdb'))


@pytest.fixture
def pdb_1hpk_roundtrip(pdb_1hpk):
    return mdt.read(pdb_1hpk.write('pdb'), format='pdb')


@pytest.mark.parametrize('mol', 'pdb_1hpk pdb_1hpk_roundtrip'.split())
def test_1hpk(request, mol):
    mol = request.getfixturevalue(mol)
    mol = mdt.interfaces.ambertools._prep_for_tleap(mol)
    for residx in (0, 21, 49, 61, 73, 78):
        residue = mol.residues[residx]
        assert residue.resname == 'CYX'
        assert residue.name[:3] == 'CYX'


@pytest.fixture
def pdb_2jaj():
    return mdt.read(get_data_path('2jaj.pdb'))


@pytest.fixture
def pdb_2jaj_roundtrip(pdb_2jaj):
    return mdt.read(pdb_2jaj.write('pdb'), format='pdb')


@pytest.mark.parametrize('mol', 'pdb_2jaj pdb_2jaj_roundtrip'.split())
def test_negative_residue_numbers_2jaj(request, mol):
    if (mol == 'pdb_2jaj_roundtrip' and
            LooseVersion(getattr(parmed, '__version__', '0.0.0')) <= LooseVersion('2.7.3')):
        pytest.xfail("This test requires ParmEd 2.7.4 (not yet released as of this writing)")

    mol = request.getfixturevalue(mol)
    res = mol.chains['B'].residues[0]
    assert res.pdbindex == -4
    assert res.index == 272
    assert res.name == 'GLY-4'


@pytest.mark.parametrize('mol', 'pdb_2jaj pdb_2jaj_roundtrip'.split())
def test_missing_residues_xtal_2jaj(request, mol):
    if mol == 'pdb_2jaj_roundtrip':
        pytest.xfail('Writing missing residue records is not yet supported.')
    mol = request.getfixturevalue(mol)
    missingres = mol.metadata.missing_residues
    for expected in MISSINGRES_2JAJ:
        assert missingres[expected[0]][expected[2]] == expected[1]


@pytest.fixture
def pdb_1pyn():
    return mdt.read(get_data_path('1pyn.pdb'))


@pytest.fixture
def pdb_1pyn_roundtrip(pdb_1pyn):
    return mdt.read(pdb_1pyn.write('pdb'), format='pdb')


@pytest.mark.parametrize('mol', 'pdb_1pyn pdb_1pyn_roundtrip'.split())
def test_numeric_residue_name_1PYN(request, mol):
    """ The ligand in this residue is named "941", which causes a little trickiness
    """
    import parmed

    mol = request.getfixturevalue(mol)
    ligand = mdt.Molecule(mol.residues[283])

    params = mdt.create_ff_parameters(ligand, charges='gasteiger')
    params._file_list['mol.lib'].put('/tmp/tmp.lib')

    contents = parmed.load_file('/tmp/tmp.lib')
    assert len(contents) == 1
    assert list(contents.keys())[0] == '941'

MISSINGRES_2JAJ = [('A', 'GLY', -4), ('A', 'PRO', -3), ('A', 'LEU', -2), ('A', 'GLY', -1),
                   ('A', 'MET', 0), ('A', 'ALA', 1), ('A', 'GLY', 2), ('A', 'LEU', 3),
                   ('A', 'GLY', 4), ('A', 'HIS', 5), ('A', 'PRO', 6), ('A', 'ALA', 7),
                   ('A', 'ALA', 32), ('A', 'LYS', 33), ('A', 'VAL', 282), ('A', 'ASP', 283),
                   ('A', 'SER', 284), ('B', 'GLY', 34), ('B', 'GLU', 35), ('B', 'ALA', 168),
                   ('B', 'ASP', 169), ('B', 'GLY', 170), ('B', 'VAL', 282), ('B', 'ASP', 283),
                   ('B', 'SER', 284)]
