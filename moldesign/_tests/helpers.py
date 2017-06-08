from builtins import range
import os

import numpy as np

import moldesign as mdt
from moldesign import units as u

DEFSTEP = 0.000005*u.angstrom


def get_data_path(f):
    """
    Returns path to a file in the ``moldesign/_tests/data/`` directory.
    """
    return os.path.join(mdt.data.PACKAGEPATH, '_tests', 'data', f)


def num_grad(mol, fn, step=DEFSTEP, atoms=None, fnargs=None, fnkwargs=None):
    """ Calculate the finite-difference gradient of a function
    """
    grad = None
    origpos = mol.positions.copy()
    if fnargs is None:
        fnargs = tuple()
    if fnkwargs is None:
        fnkwargs = dict()

    if atoms is None:
        atoms = mol.atoms

    for iatom, atom in enumerate(atoms):
        for idim in range(3):
            atom.position[idim] += step
            vplus = fn(*fnargs, **fnkwargs)
            atom.position[idim] -= 2.0 * step
            vminus = fn(*fnargs, **fnkwargs)
            mol.positions = origpos  # reset positions

            if grad is None:
                grad = np.zeros((len(atoms), 3)) * vplus.units/mol.positions.units
            grad[iatom, idim] = (vplus - vminus) / (2.0*step)

    return grad


ENERGY_TOLERANCE = 1e-8 * u.hartree

def minimization_tester(mol):
    assert 'potential_energy' not in mol.properties

    e1 = mol.calculate_potential_energy()
    p1 = mol.positions.copy()

    traj = mol.minimize()

    # check trajectory has correct initial and final positions, energies
    assert (p1 == traj.frames[0].positions).all()
    assert (mol.positions == traj.frames[-1].positions).all()
    assert abs(e1 - traj.frames[0].potential_energy) < ENERGY_TOLERANCE
    assert abs(mol.potential_energy - traj.frames[-1].potential_energy) < ENERGY_TOLERANCE

    # Force recalculation of energy to check that it's correct
    mol.calculate(use_cache=False)
    assert abs(mol.potential_energy - traj.frames[-1].potential_energy) < ENERGY_TOLERANCE


def _make_mol_with_n_hydrogens(n):
    return mdt.Molecule([mdt.Atom('H') for i in range(n)])


class ZeroEnergy(mdt.models.base.EnergyModelBase):
    """ All 0, all the time
    """
    def prep(self):
        pass

    def calculate(self):
        return dict(potential_energy=0.0*u.default.energy,
                    forces=np.zeros(self.mol.positions.shape)*u.default.force)
