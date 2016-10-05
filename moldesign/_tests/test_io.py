""" Tests for molecule creation and file i/o
"""
import collections

import pytest

import moldesign as mdt
mdt.compute.config.engine_type = 'docker'
from moldesign import units as u


@pytest.fixture
def bipyridine_sdf():
    return mdt.read('data/bipyridine.sdf')


@pytest.fixture
def bipyridine_xyz():
    return mdt.read('data/bipyridine.xyz')


@pytest.fixture
def bipyridine_mol2():
    return mdt.read('data/bipyridine.mol2')


@pytest.fixture
def bipyridine_iupac():
    return mdt.from_name('bipyridine')


@pytest.fixture
def bipyridine_smiles():
    return mdt.from_smiles('c1ccnc(c1)c2ccccn2')

ATOMDATA = {  # (symbol, valence, mass)
    1: ('H', 1, 1.008 * u.amu),
    6: ('C', 4, 12.000 * u.amu),
    7: ('N', 3, 14.003 * u.amu),
    8: ('O', 2, 15.995 * u.amu)}


@pytest.mark.parametrize('key', 'mol2 xyz sdf iupac smiles'.split())
def test_read_bipyridine_from_format(key, request):
    mol = request.getfuncargvalue('bipyridine_'+key)

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
    return mdt.read('data/ACTG.pdb')

@pytest.fixture
def dna_mmcif():
    return mdt.read('data/ACTG.cif')

@pytest.fixture
def dna_sequence():
    return mdt.build_bdna('ACTG')


@pytest.mark.parametrize('key', 'pdb mmcif sequence'.split())
def test_read_dna_from_format(key, request):
    if key == 'mmcif':
        pytest.xfail(reason='Known mmcif parser bug, fix this by 0.7.4')
    mol = request.getfuncargvalue('dna_'+key)
