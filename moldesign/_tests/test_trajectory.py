import pytest
import numpy as np

import moldesign as mdt
from moldesign import units as u

from .test_data_structures import h2_harmonic, h2_trajectory

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

    assert traj.somenumber == [1, 2, 3]
    assert traj.someletter == list('abc')
