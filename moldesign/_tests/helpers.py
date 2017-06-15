import io

from future.standard_library import install_aliases
install_aliases()
from builtins import range
from future.utils import PY2

import os
import socket

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


def _make_mol_with_n_hydrogens(n):
    return mdt.Molecule([mdt.Atom('H') for i in range(n)])


def native_str_buffer(*args, **kwargs):
    if PY2:
        return io.BytesIO(*args, **kwargs)
    else:
        return io.StringIO(*args, **kwargs)



class ZeroEnergy(mdt.models.base.EnergyModelBase):
    """ All 0, all the time
    """
    def prep(self):
        pass

    def calculate(self):
        return dict(potential_energy=0.0*u.default.energy,
                    forces=np.zeros(self.mol.positions.shape)*u.default.force)


def _internet(host="8.8.8.8", port=53, timeout=3):
    """
    Host: 8.8.8.8 (google-public-dns-a.google.com)
    OpenPort: 53/tcp
    Service: domain (DNS/TCP)

    FROM https://stackoverflow.com/a/33117579/1958900
    """
    try:
        socket.setdefaulttimeout(timeout)
        socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect((host, port))
        return True
    except Exception as ex:
        print(ex.message)
        return False

INTERNET_ON = _internet()
