import moldesign as mdt
from moldesign import units as u
import numpy as np

DEFSTEP = 0.0000005*u.angstrom


def num_grad(mol, fn, step=DEFSTEP, fnargs=None, fnkwargs=None):
    grad = None
    origpos = mol.positions.copy()
    if fnargs is None:
        fnargs = tuple()
    if fnkwargs is None:
        fnkwargs = dict()

    for iatom, atom in enumerate(mol.atoms):
        for idim in xrange(3):
            atom.position[idim] += step
            vplus = fn(*fnargs, **fnkwargs)
            atom.position[idim] -= 2.0 * step
            vminus = fn(*fnargs, **fnkwargs)
            mol.positions = origpos  # reset positions

            if grad is None:
                grad = np.zeros(mol.positions.shape) * vplus.units/mol.positions.units
            grad[iatom, idim] = (vplus - vminus) / (2.0*step)

    return grad


def _make_mol_with_n_hydrogens(n):
    return mdt.Molecule([mdt.Atom('H') for i in xrange(n)])


class ZeroEnergy(mdt.models.EnergyModelBase):
    """ All 0, all the time
    """
    def prep(self):
        pass

    def calculate(self):
        return dict(potential_energy=0.0*u.default.energy,
                    forces=np.zeros(self.mol.positions.shape)*u.default.force)