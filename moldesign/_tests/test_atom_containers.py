import itertools
import random

import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u

@pytest.fixture
def protein():
    return mdt.read('data/3AID.pdb')


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
    assert res1.distance(res2.atoms[3]) == _get_minimum_pairwise(res1, mdt.AtomList(res2.atoms[3]))


def test_center_of_mass_movement(protein):
    res = protein.residues[0]
    origpos = res.positions.copy()

    res.positions -= res.center_of_mass
    np.testing.assert_allclose(res.center_of_mass,
                               np.zeros(3))

    res.translate([1.0, 2.0, 3.0] * u.angstrom)
    np.testing.assert_allclose(res.center_of_mass,
                               [1.0, 2.0, 3.0]*u.angstrom)

    res.rotate(90 * u.degrees, [0,0,1])
    np.testing.assert_allclose(res.center_of_mass,
                               [-2.0, -1.0, 3.0]*u.angstrom)
    res.positions = origpos


def test_container_properties(protein):
    myres = random.choice(protein.residues)
    assert myres.mass == sum([atom.mass for atom in myres.atoms])
    assert myres.positions == u.array([atom.position for atom in myres.atoms])


def _get_minimum_pairwise(group1, group2):
    mindist = np.inf
    for a1, a2 in itertools.product(group1.atoms, group2.atoms):
        distance = np.sqrt(np.dot((a1.position-a2.position)**2))
        mindist = min(distance, mindist)

    return mindist
