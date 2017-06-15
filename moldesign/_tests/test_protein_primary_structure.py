""" Tests introspection and methods for protein primary structure
"""
import itertools
import pytest

import moldesign as mdt
from moldesign import units as u

from .helpers import get_data_path

fixture_types = {}


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            fixture_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@typedfixture('3AID', 'protein')
def protease_pdb():
    return mdt.read(get_data_path('3aid.pdb'))


@typedfixture('3AID', 'protein')
def protease_cif():
    return mdt.read(get_data_path('3aid.cif'))


def test_chain_types(protease_pdb):
    assert protease_pdb.chains['A'].type == 'protein'
    assert protease_pdb.chains['B'].type == 'protein'


def test_chain_iterators(protease_pdb):
    waters = list(protease_pdb.chains['A'].solvent_residues)
    waternames = [water.name for water in waters]
    assert waternames == ['HOH402', 'HOH403', 'HOH405', 'HOH406', 'HOH410']

    ligands = list(protease_pdb.chains['A'].unclassified_residues)
    assert len(ligands) == 1
    assert ligands[0].name == 'ARQ401'

    assert protease_pdb.chains['A'].get_ligand() == ligands[0]

    water_chainb = list(protease_pdb.chains['B'].solvent_residues)
    waternames_b = [water.name for water in water_chainb]
    assert waternames_b == ['HOH407', 'HOH408', 'HOH409']


def test_3aid_cif_chains(protease_cif):
    mol = protease_cif
    assert len(mol.chains) == 5
    assert mol.chains['A'].num_residues == mol.chains['B'].num_residues == 99

    assert mol.chains['C'].num_residues == 1
    assert mol.chains['C'].type == mol.chains['C'].residues[0].type == 'unknown'

    assert mol.chains['D'].type == mol.chains['D'].type == 'water'


def test_3aid_cif_separate_waters(protease_cif):
    mol = protease_cif
    assert mol.chains['D'].num_residues == 5
    assert mol.chains['E'].num_residues == 3


@pytest.mark.parametrize('fixture', fixture_types['3AID'])
def test_3aid_chain_properties(fixture, request):
    mol = request.getfixturevalue(fixture)
    for chainid in 'AB':
        c = mol.chains[chainid]
        assert c.n_terminal == c.residues['PRO1']
        assert c.c_terminal == c.residues['PHE99']
        assert c.type == 'protein'


@pytest.mark.parametrize('fixture', fixture_types['3AID'])
def test_3aid_primary_structure_access_methods(fixture, request):
    mol = request.getfixturevalue(fixture)

    a1 = mol.chains['A'].residues['GLN2'].atoms['CB']
    assert a1 is mol.atoms[a1.index]
    assert a1 is mol.chains['A'].residues.GLN2.atoms.CB
    assert a1 in mol.chains['A'].residues.GLN2
    assert 'GLN2' in mol.chains['A']
    assert hasattr(mol.chains['A'].residues.GLN2.atoms, 'CB')


@pytest.mark.parametrize('fixture', fixture_types['3AID'])
def test_3aid_atom_selection(fixture, request):
    mol = request.getfixturevalue(fixture)

    a1 = mol.chains['A'].residues['GLN2'].atoms['CB']
    a2 = mol.chains['B'].residues['LYS20'].atoms['O']
    assert abs(a1.distance(a2) - 27.206*u.angstrom) < 0.001 * u.angstrom


@pytest.mark.parametrize('fixture', fixture_types['protein'])
def test_chain_lookup_by_name_and_index(fixture, request):
    mol = request.getfixturevalue(fixture)

    for ic, chain in enumerate(mol.chains):
        assert mol.chains[chain.index] is chain
        assert mol.chains[chain.name] is chain
        assert mol.chains[ic] is chain
        # assert getattr(mol.chains, chain.name) is chain   # this functionality disabled for now, possibly forever


@pytest.mark.parametrize('fixture', fixture_types['protein'])
def test_residue_lookup_by_name_and_index(fixture, request):
    mol = request.getfixturevalue(fixture)

    for chain in mol.chains:
        for ires, residue in enumerate(chain.residues):
            assert mol.residues[residue.index] is residue
            assert chain[residue.name] is residue
            assert chain[ires] is residue
            # assert getattr(chain, residue.name) is residue  # this functionality disabled for now, possibly forever

            assert chain.residues[residue.name] is residue
            assert chain.residues[ires] is residue
            assert getattr(chain.residues, residue.name) is residue

            assert residue.chain is chain


@pytest.mark.parametrize('fixture', fixture_types['protein'])
def test_atom_lookup_by_name_and_index(fixture, request):
    mol = request.getfixturevalue(fixture)

    for residue in mol.residues:
        for iatom, atom in enumerate(residue.atoms):
            assert residue[atom.name] is atom
            assert residue[iatom] is atom
            # assert getattr(residue,atom.name) is atom  # this functionality disabled for now, possibly forever

            assert residue.atoms[atom.name] is atom
            assert residue.atoms[iatom] is atom
            assert getattr(residue.atoms, atom.name) is atom

            assert mol.atoms[atom.index] is atom

            assert atom.chain is residue.chain
            assert atom.residue is residue


@pytest.mark.parametrize('fixture', fixture_types['protein'])
def test_protein_residue_iteration(fixture, request):
    mol = request.getfixturevalue(fixture)

    assert mol.chains['A'].type == 'protein'

    firstres = mol.chains['A'].residues[0]
    for res in mol.chains['A']:
        if res.type == 'protein':
            lastres = res

    lr = firstres
    for res in mol.chains['A'][1:]:
        if res is lastres:
            break
        assert res.prev_residue is lr
        assert lr.next_residue is res
        lr = res

    lastseq = -1
    for ires, res in enumerate(mol.chains['A'].polymer_residues):
        if ires == 0:
            assert res is firstres
        assert res.pdbindex > lastseq
        lastseq = res.pdbindex
    assert res is lastres


@pytest.mark.parametrize('fixture', fixture_types['protein'])
def test_protein_residue_terminals(fixture, request):
    mol = request.getfixturevalue(fixture)

    assert mol.chains['A'].type == 'protein'

    firstres = mol.chains['A'].residues[0]
    for res in mol.chains['A']:
        if res.type == 'protein':
            lastres = res

    assert firstres.is_n_terminal
    assert lastres.is_c_terminal
    assert mol.chains['A'].n_terminal is firstres
    assert mol.chains['A'].c_terminal is lastres

    for res in mol.chains['A'][1:]:
        if res is lastres:
            break
        assert not res.is_n_terminal
        assert not res.is_c_terminal

    with pytest.raises(StopIteration):
        firstres.prev_residue

    with pytest.raises(StopIteration):
        lastres.next_residue


@pytest.mark.parametrize('fixture', fixture_types['protein'])
def test_molecule_links(fixture, request):
    mol = request.getfixturevalue(fixture)

    for obj in itertools.chain(mol.atoms, mol.residues, mol.chains):
        assert obj.molecule is mol


@pytest.mark.parametrize('fixture', fixture_types['protein'])
def test_chains_iterate_in_order(fixture, request):
    mol = request.getfixturevalue(fixture)
    _iter_index_order_tester(mol.chains)


@pytest.mark.parametrize('fixture', fixture_types['protein'])
def test_residues_iterate_in_order(fixture, request):
    mol = request.getfixturevalue(fixture)

    _iter_index_order_tester(mol.residues)

    for chain in mol.chains:
        _iter_index_order_tester(chain.residues)


@pytest.mark.parametrize('fixture', fixture_types['protein'])
def test_atoms_iterate_in_order(fixture, request):
    mol = request.getfixturevalue(fixture)

    _iter_index_order_tester(mol.atoms)

    for chain in mol.chains:
        _iter_index_order_tester(chain.atoms)

    for residue in mol.residues:
        _iter_index_order_tester(residue.atoms)


def _iter_index_order_tester(iterable):
    iterator = iter(iterable)
    lastitem = next(iterator)
    for item in iterator:
        assert item.index > lastitem.index
        lastitem = item
