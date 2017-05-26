""" Tests introspection and methods for dna primary structure
"""
import itertools
import pytest

import moldesign as mdt
from moldesign import units as u


fixture_types = {}


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            fixture_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper

@typedfixture('dna')
def bdna_actg_helix():
    dna = mdt.build_dna_helix('ACTG')
    return dna


@pytest.mark.parametrize('fixture', fixture_types['dna'])
def test_dna_chain_properties(fixture, request):
    mol = request.getfixturevalue(fixture)
    assert mol.chains['A'].type == 'dna'


@pytest.mark.parametrize('fixture', fixture_types['dna'])
def test_dna_residue_iteration(fixture, request):
    mol = request.getfixturevalue(fixture)

    assert mol.chains['A'].type == 'dna'

    firstres = mol.chains['A'].residues[0]
    for res in mol.chains['A']:
        if res.type == 'dna':
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


@pytest.mark.parametrize('fixture', fixture_types['dna'])
def test_dna_chain_terminals(fixture, request):
    mol = request.getfixturevalue(fixture)

    assert mol.chains['A'].type == 'dna'

    firstres = mol.chains['A'].residues[0]
    for res in mol.chains['A']:
        if res.type == 'dna':
            lastres = res

    assert firstres.is_5prime_end
    assert lastres.is_3prime_end
    assert mol.chains['A'].fiveprime_end is firstres
    assert mol.chains['A'].threeprime_end is lastres

    for res in mol.chains['A'][1:]:
        if res is lastres:
            break
        assert not res.is_3prime_end
        assert not res.is_5prime_end

    with pytest.raises(StopIteration):
        firstres.prev_residue

    with pytest.raises(StopIteration):
        lastres.next_residue


@pytest.mark.parametrize('fixture', fixture_types['dna'])
def test_protein_methods_on_dna_dont_work(fixture, request):
    mol = request.getfixturevalue(fixture)

    with pytest.raises(ValueError):
        mol.residues[0].is_n_terminal

    with pytest.raises(ValueError):
        mol.residues[0].is_c_terminal

    assert mol.chains[0].c_terminal is None
    assert mol.chains[0].n_terminal is None
