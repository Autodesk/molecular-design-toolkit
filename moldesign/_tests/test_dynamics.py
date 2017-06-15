import pytest

import moldesign as mdt
from moldesign import units as u

def test_shake():
    # probably there's a better way to test this, but it was already written in a notebook
    mol = mdt.from_smiles('C=C')
    monitor = mdt.DihedralMonitor(mol.atoms[0], mol.atoms[1])
    monitor.value = 90.0*u.degrees
    constraint = monitor.constrain()
    traj = mdt.Trajectory(mol)
    for i in range(50):
        oldpos = mol.positions.copy()
        monitor.atoms[0].z += 0.05*u.angstrom
        monitor.atoms[-1].x -= 0.05*u.angstrom
        mdt.geom.shake_positions(mol, oldpos)
        assert (mol.positions != oldpos).any()
        traj.new_frame()
        assert constraint.satisfied()
