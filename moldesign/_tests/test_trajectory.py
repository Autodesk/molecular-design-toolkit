import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u

from .object_fixtures import h2, h2_harmonic, h2_trajectory

def test_frames_synched_with_trajectory(h2_trajectory):
    traj = h2_trajectory

    for iframe, frame in enumerate(traj):
        assert frame.potential_energy == traj.potential_energy[iframe]
    assert iframe == traj.num_frames-1

    assert traj.potential_energy.units == traj.unit_system.energy
    assert traj.positions.units == traj.unit_system.length

    assert traj.positions.shape == (traj.num_frames, traj.mol.num_atoms, 3)
    assert traj.potential_energy.shape == (traj.num_frames,)


@pytest.fixture
def precanned_trajectory():
    a1 = mdt.Atom(6)
    a2 = mdt.Atom(1)
    a3 = mdt.Atom(7)
    mol = mdt.Molecule([a1, a2, a3])
    traj = mdt.Trajectory(mol)

    a1.x = 1.0 * u.angstrom
    a3.y = 1.0 * u.angstrom
    traj.new_frame(somenumber=1, someletter='a')

    mol.time = 1.0 * u.fs
    a1.x = 2.0 * u.angstrom
    traj.new_frame(somenumber=2, someletter='b')

    mol.time = 2.0 * u.fs
    a2.x = -1.0 * u.angstrom
    a3.x = -1.0 * u.angstrom
    traj.new_frame(somenumber=3, someletter='c')

    return traj


def test_geometry_analysis_precanned(precanned_trajectory):
    traj = precanned_trajectory
    a1, a2, a3 = traj.mol.atoms

    angles = traj.angle(a1, a2, a3)
    desired = ([90.0]*3)*u.degrees
    np.testing.assert_allclose(angles.value_in(u.degrees),
                               desired.value_in(u.degrees))

    d12 = traj.distance(a1, a2)
    desired = [1.0, 2.0, 3.0] * u.angstrom
    np.testing.assert_allclose(d12.value_in(u.angstrom),
                               desired)

    rmsd = traj.rmsd()
    desired = [0.0, 1.0/np.sqrt(3), 1.0] * u.angstrom
    np.testing.assert_allclose(rmsd.value_in(u.angstrom),
                               desired.value_in(u.angstrom))

    time = traj.time
    desired = [0.0, 1.0, 2.0] * u.fs
    np.testing.assert_allclose(time.value_in(u.fs),
                               desired.value_in(u.fs))

    assert all(traj.somenumber == np.array([1, 2, 3]))
    assert traj.someletter == list('abc')


def test_frame_to_molecule_conversion(precanned_trajectory):
    traj = precanned_trajectory

    f0_position = u.angstrom*[[1, 0, 0],
                              [0, 0, 0],
                              [0, 1, 0]]

    f2_position = u.angstrom*[[2, 0, 0],
                              [-1, 0, 0],
                              [-1, 1, 0]]

    f0 = traj.frames[0]
    mol = f0.as_molecule()
    assert mol.same_topology(traj.mol)
    assert (mol.positions == f0_position).all()
    assert mol.time == 0.0 * u.fs

    # test ability to directly write a frame
    readmol = mdt.read(f0.write('xyz'), format='xyz')
    np.testing.assert_allclose(readmol.positions.value_in(u.angstrom),
                               f0_position.value_in(u.angstrom))


    m2 = traj.frames[-1].as_molecule()
    assert m2.same_topology(mol)
    assert (m2.positions == f2_position).all()
    assert m2.time == 2.0 * u.fs


def test_property_backfill(precanned_trajectory):
    traj = precanned_trajectory
    oldnumframes = len(traj)

    traj.mol.time += 1.0*u.fs
    traj.new_frame(somenewthing=5)

    assert traj.somenewthing == [None] * oldnumframes + [5]