from builtins import range
import itertools
import random

import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u

from . import helpers


registered_types = {}


def typedfixture(*types, **kwargs):
    """This is a decorator that lets us associate fixtures with one or more arbitrary types.
    We'll later use this type to determine what tests to run on the result"""

    def fixture_wrapper(func):
        for t in types:
            registered_types.setdefault(t, []).append(func.__name__)
        return pytest.fixture(**kwargs)(func)

    return fixture_wrapper


@typedfixture('container')
def protein():
    return mdt.read(helpers.get_data_path('3aid.pdb'))


@pytest.fixture
def atom(protein):
    return random.choice(protein.atoms)


@typedfixture('container')
def chain(protein):
    return random.choice(protein.chains)


@typedfixture('container')
def residue(protein):
    return random.choice(protein.residues)


@typedfixture('container')
def atomlist(protein):
    return mdt.AtomList(random.sample(protein.atoms, 10))


@typedfixture('container')
def small_molecule():
    return mdt.from_smiles('c1ccccc1')


def test_self_distances_are_0(protein):
    res = protein.residues[0]

    assert protein.distance(protein.atoms[0]) == 0.0
    assert protein.residues[0].distance(protein) == 0.0
    assert protein.atoms[0].distance(protein.atoms[0]) == 0.0
    assert res.distance(res) == 0.0
    assert res.distance(res.chain) == 0.0


def test_distance_is_minimum_pairwise(protein):
    res1 = protein.residues[0]
    res2 = protein.residues[2]

    assert res1.distance(res2) == _get_minimum_pairwise(res1, res2)
    assert res1.distance(res2.atoms[3]) == _get_minimum_pairwise(res1,
                                                                 mdt.AtomList([res2.atoms[3]]))


@pytest.mark.parametrize(['f1', 'f2'],
                         itertools.product(registered_types['container'],
                                           registered_types['container']))
def test_pairwise_distance_arrays(f1, f2, request):
    o1 = request.getfixturevalue(f1)
    o2 = request.getfixturevalue(f2)

    array = o1.calc_distance_array(o2)

    if o1.num_atoms * o2.num_atoms > 250:  # stochastically test larger matrices
        pairs = ((random.randrange(0, o1.num_atoms), random.randrange(0, o2.num_atoms))
                 for i in range(250))
    else:
        pairs = itertools.product(range(o1.num_atoms), range(o2.num_atoms))

    for iatom, jatom in pairs:
        assert o1.atoms[iatom].distance(o2.atoms[jatom]) == array[iatom, jatom]


@pytest.mark.parametrize('fixturename', registered_types['container'])
def test_center_of_mass_movement(fixturename, request):
    obj = request.getfixturevalue(fixturename)
    origpos = obj.positions.copy()

    obj.positions -= obj.center_of_mass
    np.testing.assert_allclose(obj.center_of_mass,
                               np.zeros(3),
                               atol=1e-10)

    obj.translate([1.0, 2.0, 3.0] * u.angstrom)
    np.testing.assert_allclose(obj.center_of_mass,
                               [1.0, 2.0, 3.0]*u.angstrom)

    obj.rotate(90 * u.degrees, [0,0,1])
    np.testing.assert_allclose(obj.center_of_mass,
                               [1.0, 2.0, 3.0]*u.angstrom)

    obj.rotate(90 * u.degrees, [0,0,1], center=[0,0,0]*u.angstrom)
    np.testing.assert_allclose(obj.center_of_mass,
                               [-2.0, 1.0, 3.0]*u.angstrom)
    obj.positions = origpos


@pytest.mark.parametrize('fixturename', registered_types['container'])
def test_container_properties(fixturename, request):
    obj = request.getfixturevalue(fixturename)
    assert obj.mass == sum([atom.mass for atom in obj.atoms])
    np.testing.assert_array_equal(obj.positions.defunits(),
                                  u.array([atom.position for atom in obj.atoms]).defunits())
    assert obj.num_atoms == len(obj.atoms)


@pytest.mark.parametrize('fixturename', registered_types['container'])
def test_position_links(fixturename, request):
    obj = request.getfixturevalue(fixturename)

    np.testing.assert_array_equal(obj.positions[0, :],
                                  obj.atoms[0].position)
    obj.positions[0, :] *= 2.0
    np.testing.assert_array_equal(obj.positions[0, :],
                                  obj.atoms[0].position)


def _get_minimum_pairwise(group1, group2):
    mindist = np.inf * u.angstrom
    for a1, a2 in itertools.product(group1.atoms, group2.atoms):
        distance = (a1.position-a2.position).norm()
        mindist = min(distance, mindist)

    return mindist


@pytest.mark.parametrize('fixturename', ['atom', 'residue', 'atomlist', 'small_molecule'])
def test_atoms_within(fixturename, request):
    obj = request.getfixturevalue(fixturename)

    if fixturename == 'atom':
        myatoms = {obj}
        mol = obj.molecule
    else:
        myatoms = set(obj.atoms)
        mol = obj.atoms[0].molecule

    assert len(obj.atoms_within(0.0*u.angstrom)) == 0

    within5 = set(obj.atoms_within(5.0*u.angstrom))
    within5_self = set(obj.atoms_within(5.0*u.angstrom, include_self=True))

    assert myatoms.issubset(within5_self)
    assert within5.union(myatoms) == within5_self

    for atom in mol.atoms:
        if atom in myatoms:
            assert atom not in within5
        elif atom.distance(obj) <= 5.0*u.angstrom:
            assert atom in within5
        else:
            assert atom not in within5


@pytest.mark.parametrize('fixturename', ['atom', 'residue', 'atomlist'])
def test_residues_within(fixturename, request):
    obj = request.getfixturevalue(fixturename)

    if fixturename == 'atom':
        mol = obj.molecule
    else:
        mol = obj.atoms[0].molecule

    assert len(obj.residues_within(0.0*u.angstrom)) == 0

    within5 = set(obj.residues_within(5.0*u.angstrom))
    within5_self = set(obj.residues_within(5.0*u.angstrom, include_self=True))

    if fixturename == 'residue':
        assert within5.union([obj]) == within5_self
    else:
        assert within5 == within5_self

    for residue in mol.residues:
        if residue == obj:
            assert residue in within5_self
            assert residue not in within5
        elif residue.distance(obj) <= 5.0*u.angstrom:
            assert residue in within5
        else:
            assert residue not in within5


def test_setlike_atomlist_methods(protein):
    l1 = protein.atoms[:10]
    assert isinstance(l1, mdt.AtomList)
    l2 = protein.atoms[5:15]

    assert l1.union(l2) == protein.atoms[:15]
    assert l2.union(l1) == protein.atoms[5:15] + protein.atoms[:5]

    interx = l1.intersection(l2)
    assert interx == protein.atoms[5:10]
    assert l2.intersection(l1) == interx

    assert l1 - l2 == protein.atoms[:5]
    assert l2 - l1 == protein.atoms[10:15]

    assert (l1 + l2).unique() == protein.atoms[:15]


def test_get_atoms(protein):
    protein_atoms = set(protein.get_atoms('protein'))
    water_atoms = set(protein.get_atoms('water'))
    carbon_atoms = set(protein.get_atoms(symbol='C'))
    gly_alpha_carbons = set(protein.get_atoms(resname='GLY', name='CA'))
    for atom in protein.atoms:
        assert (atom.residue.type == 'protein') == (atom in protein_atoms)
        assert (atom.residue.resname == 'HOH') == (atom in water_atoms)
        assert (atom.atnum == 6) == (atom in carbon_atoms)
        assert (atom.name == 'CA' and atom.residue.resname == 'GLY') == (
            atom in gly_alpha_carbons)

