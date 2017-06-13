import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u
from moldesign.external import transformations


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
                'C2': 1}
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


@pytest.mark.parametrize('fixturename', 'linear linear_inversion planar'.split())
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
        assert_isomorphic(mol, positions, transformed, tolerance=elem.csm)

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
            assert abs(angle.value_in(u.radian) - theta) <= 1e-6

def assert_isomorphic(mol, transformed, positions, tolerance):
    if tolerance == 0.0:
        tolerance = 1.0e-8 * u.angstrom

    tol2 = (tolerance * 1.1) ** 2

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