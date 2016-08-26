""" Tests introspection and methods for protein primary structure
"""

import pytest

import moldesign as mdt
from moldesign import units as u


registered_types = {}


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@typedfixture('3AID', 'protein')
def protease_pdb():
    return mdt.read('data/3aid.pdb')


@typedfixture('3AID', 'protein')
def protease_cif():
    return mdt.read('data/3aid.cif')


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


@pytest.mark.parametrize('objkey', registered_types['3AID'])
def test_3aid_chain_properties(objkey, request):
    mol = request.getfuncargvalue(objkey)
    for chainid in 'AB':
        c = mol.chains[chainid]
        assert c.n_terminal == c.residues['PRO1']
        assert c.c_terminal == c.residues['PHE99']
        assert c.type == 'protein'


@pytest.mark.parametrize('objkey', registered_types['3AID'])
def test_3aid_atom_selection(objkey, request):
    mol = request.getfuncargvalue(objkey)

    a1 = mol.chains['A'].residues['GLN2'].atoms['CB']
    a2 = mol.chains['B'].residues['LYS20'].atoms['O']
    assert abs(a1.distance(a2) - 27.206*u.angstrom) < 0.001 * u.angstrom


@pytest.mark.parametrize('objkey', registered_types['protein'])
def test_chain_lookup_by_name_and_index(objkey, request):
    mol = request.getfuncargvalue(objkey)

    for chain in mol.chains:
        assert mol.chains[chain.index] is chain
        assert mol.chains[chain.name] is chain


@pytest.mark.parametrize('objkey', registered_types['protein'])
def test_residue_lookup_by_name_and_index(objkey, request):
    mol = request.getfuncargvalue(objkey)

    for chain in mol.chains:
        for residue in chain.residues:
            assert mol.residues[residue.index] is residue
            assert chain[residue.name] is residue

            assert residue.chain is chain


@pytest.mark.parametrize('objkey', registered_types['protein'])
def test_atom_lookup_by_name_and_index(objkey, request):
    mol = request.getfuncargvalue(objkey)

    for residue in mol.residues:
        for atom in residue.atoms:
            assert residue[atom.name] is atom
            assert mol.atoms[atom.index] is atom

            assert atom.chain is residue.chain
            assert atom.residue is residue


@pytest.mark.parametrize('objkey', registered_types['protein'])
def test_chains_iterate_in_order(objkey, request):
    mol = request.getfuncargvalue(objkey)
    _iter_index_order_tester(mol.chains)


@pytest.mark.parametrize('objkey', registered_types['protein'])
def test_residues_iterate_in_order(objkey, request):
    mol = request.getfuncargvalue(objkey)
    _iter_index_order_tester(mol.residues)

    for chain in mol.chains:
        _iter_index_order_tester(chain.residues)


@pytest.mark.parametrize('objkey', registered_types['protein'])
def test_atoms_iterate_in_order(objkey, request):
    mol = request.getfuncargvalue(objkey)
    _iter_index_order_tester(mol.atoms)

    for chain in mol.chains:
        _iter_index_order_tester(chain.atoms)

    for residue in mol.residues:
        _iter_index_order_tester(residue.atoms)


def _iter_index_order_tester(iterable):
    iterator = iter(iterable)
    lastitem = iterator.next()
    for item in iterator:
        assert item.index > lastitem.index
        lastitem = item
