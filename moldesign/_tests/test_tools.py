""" Tests topology manipulation tools
"""
import collections

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


@pytest.fixture
def ammonium_nocharge():
    return mdt.from_smiles('[NH4]')


@pytest.fixture
def ammonium_charged():
    return mdt.from_smiles('[NH4+]')


@pytest.mark.parametrize('objkey',
                         ['ammonium_nocharge', 'ammonium_charged'])
def test_ammonium_formal_charge(objkey, request):
    mol = request.getfuncargvalue(objkey)
    mdt.assign_formal_charges(mol)
    assert mol.charge == 1 * u.q_e

    for atom in mol.atoms:
        if atom.atnum == 7:
            assert atom.formal_charge == 1 * u.q_e
        else:
            assert atom.atnum == 1
            assert atom.formal_charge == 0 * u.q_e


@pytest.fixture
def c2_no_hydrogen_from_smiles():
    mymol = mdt.from_smiles('[CH0][CH0]')
    return mymol


def test_c2_no_hydrogen_from_smiles(c2_no_hydrogen_from_smiles):
    mymol = c2_no_hydrogen_from_smiles

    atomcounts = collections.Counter(atom.element for atom in mymol.atoms)
    assert atomcounts['C'] == 2
    assert len(atomcounts) == 1
    assert mymol.num_bonds == 1
    assert mymol.num_atoms == 2
    bonds = list(mymol.bonds)
    assert len(bonds) == 1
    b = bonds[0]
    assert b.order == 1
    assert b.a1.index == 0
    assert b.a2.index == 1


def test_add_hydrogen_to_c2(c2_no_hydrogen_from_smiles):
    newmol = mdt.add_hydrogen(c2_no_hydrogen_from_smiles)
    atomcounts = collections.Counter(atom.element for atom in newmol.atoms)
    assert newmol.num_atoms == 8
    assert atomcounts['C'] == 2
    assert atomcounts['H'] == 6
    assert len(atomcounts) == 2
    assert newmol.num_bonds == 7
    for atom, bondgraph in newmol.bond_graph.iteritems():
        if atom.atnum == 1:
            assert len(bondgraph) == 1
            assert bondgraph.keys()[0].elem == 'C'
            assert bondgraph.values()[0] == 1
        else:
            assert atom.atnum == 6
            assert len(bondgraph) == 4
            for nbr in bondgraph:
                assert bondgraph[nbr] == 1




