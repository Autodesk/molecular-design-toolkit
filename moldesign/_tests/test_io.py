""" Tests for molecule creation and file i/o
"""
from builtins import str
import collections

import numpy
import pytest

import moldesign as mdt
mdt.compute.config.engine_type = 'docker'
from moldesign import units as u

from .helpers import get_data_path


@pytest.fixture
def bipyridine_sdf():
    return mdt.read(get_data_path('bipyridine.sdf'))


@pytest.fixture
def bipyridine_xyz():
    return mdt.read(get_data_path('bipyridine.xyz'))


@pytest.fixture
def bipyridine_mol2():
    return mdt.read(get_data_path('bipyridine.mol2'))


@pytest.fixture
def bipyridine_iupac():
    return mdt.from_name('bipyridine')

@pytest.fixture
def bipyridine_inchi():
    return mdt.from_inchi('InChI=1S/C10H8N2/c1-3-7-11-9(5-1)10-6-2-4-8-12-10/h1-8H')


@pytest.fixture
def bipyridine_smiles():
    return mdt.from_smiles('c1ccnc(c1)c2ccccn2')

ATOMDATA = {  # (symbol, valence, mass)
    1: ('H', 1, 1.008 * u.amu),
    6: ('C', 4, 12.000 * u.amu),
    7: ('N', 3, 14.003 * u.amu),
    8: ('O', 2, 15.995 * u.amu)}


@pytest.mark.parametrize('key', 'iupac smiles inchi xyz sdf'.split())
def test_auto_unique_atom_names(key, request):
    mol = request.getfixturevalue('bipyridine_'+key)

    atomnames = set(atom.name for atom in mol.atoms)
    assert len(atomnames) == mol.num_atoms


def test_atom_names_preserved_from_input_file_mol2(bipyridine_mol2):
    mol = bipyridine_mol2
    for atom in mol.atoms:
        assert atom.name == atom.symbol + str(atom.index)


@pytest.mark.parametrize('key', 'mol2 xyz sdf iupac smiles inchi'.split())
def test_read_bipyridine_from_format(key, request):
    mol = request.getfixturevalue('bipyridine_'+key)

    atomcounts = collections.Counter(atom.symbol for atom in mol.atoms)
    assert len(atomcounts) == 3
    assert atomcounts['C'] == 10
    assert atomcounts['N'] == 2
    assert atomcounts['H'] == 8

    assert mol.charge == 0
    assert abs(mol.mass - 156.069*u.amu) < 0.001 * u.amu
    for atom in mol.atoms:
        assert atom.formal_charge == 0.0
        symb, val, mss = ATOMDATA[atom.atnum]
        assert atom.symbol == symb
        assert atom.valence == val
        assert abs(atom.mass - mss) < 0.001 * u.amu

    assert mol.num_bonds == 21
    bondorders = collections.Counter(bond.order for bond in mol.bonds)
    assert bondorders[2] == 6
    assert bondorders[1] == 15
    assert len(bondorders) == 2


@pytest.fixture
def dna_pdb():
    return mdt.read(get_data_path('ACTG.pdb'))


@pytest.fixture
def dna_mmcif():
    return mdt.read(get_data_path('ACTG.cif'))


@pytest.fixture
def dna_sequence():
    return mdt.build_bdna('ACTG')


@pytest.fixture
def pdb_1kbu():
    return mdt.read(get_data_path('1KBU.pdb'))


@pytest.fixture
def mmcif_1kbu():
    return mdt.read(get_data_path('1KBU.cif'))


@pytest.mark.parametrize('key', 'pdb mmcif sequence'.split())
def test_read_dna_from_format(key, request):
    if key == 'mmcif':
        pytest.xfail(reason='Known mmcif parser bug, fix this by 0.7.4')
    mol = request.getfixturevalue('dna_'+key)


@pytest.mark.parametrize('key', 'mmcif pdb'.split())
def test_1kbu_assembly_data(key, request):
    mol = request.getfixturevalue('%s_1kbu' % key)

    assert len(mol.properties.bioassemblies) == 1
    assert '1' in mol.properties.bioassemblies
    assembly = mol.properties.bioassemblies['1']

    assert len(assembly.transforms) == 2
    assert set(assembly.chains) == set(c.name for c in mol.chains)

    # first transform is identity
    numpy.testing.assert_allclose(assembly.transforms[0],
                                  numpy.identity(4))

    # second transform's rotation is unitary
    rot = assembly.transforms[1][:3,:3]
    numpy.testing.assert_allclose(rot.dot(rot.T),
                                  numpy.identity(3))


@pytest.mark.parametrize('key', 'mmcif pdb'.split())
def test_1kbu_assembly_build(key, request):
    asym = request.getfixturevalue('%s_1kbu' % key)

    original = mdt.Molecule(asym)

    assembly = asym.properties.bioassemblies['1']

    rot = assembly.transforms[1][:3,:3]
    move = assembly.transforms[1][:3,3] * u.angstrom

    mol = mdt.build_assembly(asym, 1)
    assert mol.num_chains == 2 * asym.num_chains

    # test that original is unaffected
    assert original.is_identical(asym)

    testchain = assembly.chains[0]
    new_chain_pos = mol.chains[testchain].positions.T.ldot(rot).T + move[None, :]
    numpy.testing.assert_allclose(new_chain_pos.defunits_value(),
                                  mol.chains[asym.num_chains].positions.defunits_value())


@pytest.mark.parametrize('fmt', 'smiles pdb mol2 sdf inchi mmcif pkl'.split())
def test_topology_preserved_in_serialization(bipyridine_smiles, fmt):
    """ Test that bond topology is preserved even if it doesn't make sense from distances
    """
    if fmt != 'pkl':
        pytest.xfail("We are currently unable to get an unambiguous representation of a molecular "
                     "sructure with ANY current file formats or parsers.")
    mol = bipyridine_smiles.copy()  # don't screw up the fixture object
    mol.bond_graph[mol.atoms[3]][mol.atoms[5]] = 3
    mol.bond_graph[mol.atoms[5]][mol.atoms[3]] = 3
    mol.atoms[3].x += 10.0 * u.angstrom

    newmol = mdt.read(mol.write(format=fmt), format=fmt)
    assert mol.same_bonds(newmol, verbose=True)
