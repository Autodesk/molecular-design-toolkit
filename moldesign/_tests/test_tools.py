""" Tests topology manipulation tools
"""
import collections
import random

import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u

from .helpers import get_data_path


registered_types = {}


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@pytest.fixture(scope='function')
def linear_mol_and_pmi():
    mol = mdt.from_smiles('C#CC#CC#C')
    pmi = mdt.tools.PrincipalMomentsOfInertia(mol)
    return mol, pmi


def test_pmi_orthonormal(linear_mol_and_pmi):
    pmi, mol = linear_mol_and_pmi
    for ivec in range(3):
        assert abs(1.0 - mdt.mathutils.norm(pmi.evecs[ivec])) < 1e-12
        for jvec in range(ivec, 3):
            assert abs(pmi.evecs[ivec].dot(pmi.evecs[jvec])) < 1e-12


def test_principal_moment_of_inertia_reorientation(linear_mol_and_pmi):
    # note that this will fail if the generated polymer is not perfectly linear
    mol, pmi = linear_mol_and_pmi

    original_distmat = mol.calc_distance_array().defunits_value()

    for deg in np.linspace(10, 120, 6) * u.degree:
        # Set up a new random reorientation every time
        axis = [random.uniform(-1.0, 1.0),
                random.uniform(-1.0, 1.0),
                random.uniform(0.4, 1.0)]  # make sure to rotate off of x-axis
        translation = 10.0 * (np.random.rand(3)-0.5) * u.angstrom
        mol.translate(translation)
        mol.rotate(axis=axis, angle=deg)

        # calculate PMI and reorient
        pmi = mdt.tools.PrincipalMomentsOfInertia(mol)
        pmi.reorient(mol)

        # Check that everything lies along the x-axis
        np.testing.assert_allclose(mol.positions[:,1:], 0.0)

        # Check that distance matrix was unchanged under the transformation
        newdistmat = mol.calc_distance_array()
        np.testing.assert_allclose(newdistmat.defunits_value,
                                   original_distmat)


@pytest.fixture
def ammonium_nocharge():
    return mdt.from_smiles('[NH4]')


@pytest.fixture
def ammonium_charged():
    return mdt.from_smiles('[NH4+]')


@pytest.mark.parametrize('objkey',
                         ['ammonium_nocharge', 'ammonium_charged'])
def test_ammonium_formal_charge(objkey, request):
    mol = request.getfixturevalue(objkey)
    mdt.assign_formal_charges(mol)
    assert mol.charge == 1 * u.q_e

    for atom in mol.atoms:
        if atom.atnum == 7:
            assert atom.formal_charge == 1 * u.q_e
        else:
            assert atom.atnum == 1
            assert atom.formal_charge == 0 * u.q_e


def test_set_hybridization_and_saturate():
    # Creates just the carbons of ethylene, expects the routine to figure out the rest
    atom1 = mdt.Atom(6)
    atom2 = mdt.Atom(6)
    atom2.x = 1.35 * u.angstrom
    atom1.bond_to(atom2, 1)
    mol = mdt.Molecule([atom1, atom2])
    newmol = mdt.set_hybridization_and_saturate(mol)
    pytest.xfail('This is apparently broken')
    assert newmol.num_atoms == 6
    assert newmol.atoms[0].bond_graph[atom1] == 2
    assert len(newmol.get_atoms(atnum=1)) == 4


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
    for atom, bondgraph in newmol.bond_graph.items():
        if atom.atnum == 1:
            assert len(bondgraph) == 1
            assert list(bondgraph.keys())[0].elem == 'C'
            assert list(bondgraph.values())[0] == 1
        else:
            assert atom.atnum == 6
            assert len(bondgraph) == 4
            for nbr in bondgraph:
                assert bondgraph[nbr] == 1


def test_split_chains_3p3k():
    mol = mdt.read(get_data_path('3p3k.pdb'))
    assert mol.num_chains == 1  # not really testing this, just for sanity's sake

    newmol = mdt.split_chains(mol)
    assert newmol.num_chains == 3

    assert newmol.chains['A'] is newmol.chains[0]
    assert newmol.chains['A'].type == 'protein'

    assert newmol.chains['B'] is newmol.chains[1]
    assert newmol.chains['B'].type == 'protein'

    assert newmol.chains['C'] is newmol.chains[2]
    assert newmol.chains['C'].type == 'water'

