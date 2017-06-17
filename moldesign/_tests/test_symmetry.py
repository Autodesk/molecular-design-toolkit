import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u
from moldesign.external import transformations

from .molecule_fixtures import benzene


@pytest.fixture
def linear():
    expected = {'C1': 1,
                'Cinf_v': 1}

    a1 = mdt.Atom(1)
    a1.position = [-1, 0, 0] * u.angstrom
    a2 = mdt.Atom(2)
    a2.position = [1, 0, 0] * u.angstrom

    return mdt.Molecule([a1, a2], name='linear_h2'), expected, True


@pytest.fixture
def linear_inversion(linear):
    mol = linear[0].copy()

    expected = {'C1': 1,
                'Cinf_v': 1,
                'Cs': 1,
                'Ci': 1}
    mol.atoms[0].atnum = mol.atoms[1].atnum
    mol.atoms[0].mass = mol.atoms[1].mass

    return mol, expected, True


@pytest.fixture
def planar():
    expected = {'C1': 1,
                'Cs': 1}

    a1 = mdt.Atom(1)
    a1.position = [-1, 0, 0] * u.angstrom
    a2 = mdt.Atom(2)
    a2.position = [0.9, 0.2, 0] * u.angstrom
    a3 = mdt.Atom(3)
    a3.position = [-0.9, 3.2, 0] * u.angstrom

    return mdt.Molecule([a1, a2, a3], name='planar_h3'), expected, True


@pytest.fixture
def almost_planar(planar):
    atoms = planar[0].atoms
    atoms.append(mdt.Atom(4))
    atoms[-1].position = [0, 0, 0.005] * u.angstrom
    return mdt.Molecule(atoms, name='almost_planar')


def test_approximate_symmetry(almost_planar):
    mol = almost_planar.copy()

    symmetries = mdt.get_symmetry(almost_planar)

    assert symmetries.rms > 0
    assert len(symmetries.exact) == 1
    assert symmetries.exact[0].symbol == 'C1'

    assert len(symmetries.approximate) == 1
    assert symmetries.approximate[0].symbol == 'Cs'
    cs = symmetries.groups['Cs'][0]

    mol.positions = symmetries.get_symmetrized_coords(cs)
    newsymm = mdt.get_symmetry(mol)
    assert len(newsymm.exact) == 2
    assert len(newsymm.approximate) == 0


@pytest.fixture
def benz(benzene):
    # Note: we're counting on the geometry generation algorithm to give proper symmetries
    expected = {'C1': 1,
                'C6': 1,
                'C2': 7,
                'C3': 1,
                'Ci': 1,
                'Cs': 7,
                'S6': 1,
                'S3': 1}
    return benzene, expected, False


@pytest.mark.parametrize('fixturename', 'linear linear_inversion planar benz'.split())
def test_expected_symmetries(fixturename, request):
    mol, expected, allexact = request.getfixturevalue(fixturename)
    symmetries = mdt.get_symmetry(mol)

    if allexact:
        assert abs(symmetries.rms) <= 1e-10 * u.angstrom

    found = {}
    for elem in symmetries.elems:
        if elem.symbol not in found:
            found[elem.symbol] = 0
        found[elem.symbol] += 1

        if allexact:
            assert abs(elem.csm) <= 1e-10 * u.angstrom
            assert abs(elem.max_diff) <= 1e-10 * u.angstrom

        positions = symmetries.orientation.copy()
        transformed = positions.dot(elem.matrix)
        assert_isomorphic(mol, positions, transformed, elem.max_diff)

        # check that transformation matrix is unitary
        assert_identity_matrix(elem.matrix.dot(elem.matrix.T))

        # check matrix transform properties
        if elem.symbol == 'C1':
            assert_identity_matrix(elem.matrix)
        elif elem.symbol == 'Cs':
            assert_reflection_matrix(elem.matrix)
        elif elem.symbol == 'Ci':
            assert_identity_matrix(-1 * elem.matrix)
        elif len(elem.symbol) == 2 and elem.symbol[0] == 'C' and elem.symbol[1].isdigit():
            nfold = int(elem.symbol[1])
            assert_rotation_matrix(elem.matrix, (360.0*u.degrees)/nfold)
    assert expected == found


def assert_identity_matrix(mat):
    np.testing.assert_allclose(mat, np.identity(3), atol=1e-10)


def assert_reflection_matrix(mat):
    mat4 = _asmat4(mat)
    try:
        transformations.reflection_from_matrix(mat4)
    except ValueError:
        assert False, 'not a reflection matrix'


def assert_rotation_matrix(mat, angle=None):
    mat4 = _asmat4(mat)
    try:
        theta, direc, point = transformations.rotation_from_matrix(mat4)
    except ValueError:
        assert False, 'not a rotation matrix'

    else:
        assert (point[:3] == 0).all()  # would be WEIRD if this fails

        if angle is not None:
            assert abs(angle.value_in(u.radian) - abs(theta)) <= 1e-6


def assert_isomorphic(mol, transformed, positions, maxdiff):
    if maxdiff == 0.0:
        tolerance = 1.0e-4 * u.angstrom
    else:
        tolerance = maxdiff * 4

    tol2 = tolerance ** 2

    available = {atom: pos for atom, pos in zip(mol.atoms, positions)}

    for atom, newpos in zip(mol.atoms, transformed):
        closest_atom = None
        closest = np.inf * u.angstrom**2
        for availatom, pos in available.items():
            dist2 = np.sum((pos-newpos)**2)
            if dist2 < closest:
                closest_atom = availatom
                closest = dist2

        assert atom.atnum == closest_atom.atnum
        assert atom.mass == closest_atom.mass
        assert closest <= tol2


def _asmat4(mat):
    mat4 = np.identity(4)
    mat4[:3,:3] = mat
    return mat4