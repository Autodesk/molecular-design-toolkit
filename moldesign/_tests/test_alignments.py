import random

import numpy as np
import pytest

import moldesign as mdt
from moldesign import geom
from moldesign.mathutils import normalized
from moldesign import units as u

from .molecule_fixtures import *
from .helpers import assert_almost_equal


__PYTEST_MARK__ = 'internal'  # mark all tests in this module with this label (see ./conftest.py)


def test_bond_alignment_on_axis(benzene):
    mol = benzene.copy()
    directions = ['x', 'y', 'z',
                  [1,2,3.0],
                  [0,1,0],
                  [0.1, 0.1, 0.1] * u.angstrom]

    for i, dir in enumerate(directions):
        bond = mdt.Bond(*random.sample(mol.atoms, 2))
        center = (i % 2) == 0.0
        bond.align(dir, centered=center)

        if center:
            np.testing.assert_allclose(bond.midpoint, np.zeros(3),
                                       atol=1.e-12)
            np.testing.assert_allclose(bond.a1.position.defunits_value(),
                                       -bond.a2.position.defunits_value(),
                                       atol=1e-10)
        if isinstance(dir, str):
            if dir == 'x':
                d = np.array([1.0, 0, 0])
            elif dir == 'y':
                d = np.array([0.0, 1.0, 0.0])
            elif dir == 'z':
                d = np.array([0.0, 0.0, 1.0])
            else:
                raise WtfError()
        else:
            d = normalized(u.array(dir))

        newvec = (bond.a2.position - bond.a1.position).normalized()
        assert abs(1.0 - d.dot(newvec)) < 1e-10


def test_align_two_bonds(benzene, h2):
    h2 = h2.copy()
    alignbond = list(h2.bonds)[0]
    one_was_different = False

    for targetbond in benzene.bonds:

        # sanity checks before we start:
        assert not np.allclose(alignbond.midpoint.defunits_value(),
                               targetbond.midpoint.defunits_value())
        vec1, vec2 = (normalized(b.a2.position - b.a1.position)
                      for b in (alignbond, targetbond))
        one_was_different = (one_was_different or vec1.dot(vec2) < 0.8)

        alignbond.align(targetbond)

        np.testing.assert_allclose(alignbond.midpoint, targetbond.midpoint)
        vec1, vec2 = (normalized(b.a2.position - b.a1.position)
                      for b in (alignbond, targetbond))
        assert (1.0 - vec1.dot(vec2)) < 1e-8


@pytest.mark.parametrize('fixturename', 'ligand3aid pdb3aid benzene small_molecule pdb1yu8'.split())
def test_pmi_is_orthonormal(request, fixturename):
    # test a bunch of molecules
    mol = request.getfixturevalue(fixturename).copy()
    pmi = geom.alignment.PrincipalMomentsOfInertia(mol)
    for ivec in range(3):
        assert abs(1.0 - mdt.mathutils.norm(pmi.evecs[ivec])) < 1e-12
        for jvec in range(ivec+1, 3):
            assert abs(pmi.evecs[ivec].dot(pmi.evecs[jvec])) < 1e-12


@pytest.mark.parametrize('fixturename', 'ligand3aid pdb3aid benzene small_molecule pdb1yu8'.split())
def test_pmi_translational_rotational_invariance(request, fixturename):
    # currently only test the eigenvalues
    mol = request.getfixturevalue(fixturename).copy()
    pmi = geom.alignment.PrincipalMomentsOfInertia(mol)

    _randomize_orientation(mol)
    pmi2 = geom.alignment.PrincipalMomentsOfInertia(mol)
    assert_almost_equal(pmi.moments, pmi2.moments)


def _randomize_orientation(mol):
    translation = 10.0*(np.random.rand(3)-0.5)*u.angstrom
    mol.translate(translation)

    axis = [random.uniform(-1.0, 1.0),
            random.uniform(-1.0, 1.0),
            random.uniform(0.4, 1.0)]  # make sure to rotate off of x-axis
    mol.rotate(axis=axis, angle=random.randrange(15.0, 170.0) * u.degrees)


@pytest.mark.screening
def test_pmi_reorientation_on_linear_molecule():
    # note that this will fail if the generated polymer is not perfectly linear
    mol = mdt.from_smiles('C#CC#CC#C')
    original_distmat = mol.calc_distance_array().defunits_value()

    for i in range(5):
        _randomize_orientation(mol)

        # calculate PMI and reorient
        pmi = geom.alignment.PrincipalMomentsOfInertia(mol)
        pmi.reorient(mol)

        # Check that everything lies along the z-axis
        np.testing.assert_allclose(mol.positions[:,:2], 0.0, atol=1e-12)

        # Check that distance matrix was unchanged under the transformation
        newdistmat = mol.calc_distance_array()
        np.testing.assert_allclose(newdistmat.defunits_value(),
                                   original_distmat)


def test_pmi_reorientation_on_benzene(benzene):
    # note that this will fail if the generated polymer is not perfectly linear
    mol = benzene.copy()
    original_distmat = mol.calc_distance_array().defunits_value()

    for i in range(5):
        _randomize_orientation(mol)

        # calculate PMI and reorient
        pmi = geom.alignment.PrincipalMomentsOfInertia(mol)
        pmi.reorient(mol)

        # Check that molecule is _roughly_ in x-y plane
        np.testing.assert_allclose(mol.positions[:,0].value_in(u.angstrom), 0.0,
                                   atol=0.1)

        # Check that distance matrix was unchanged under the transformation
        newdistmat = mol.calc_distance_array()
        np.testing.assert_allclose(newdistmat.defunits_value(),
                                   original_distmat)

