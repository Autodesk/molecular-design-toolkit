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


ENERGY_TOLERANCE = 1e-8 * u.hartree


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
        try:
            socket.socket(socket.AF_INET, socket.SOCK_STREAM).connect((host, port))
        except OSError:
            return False
        return True
    except Exception as ex:
        print(ex.message)
        return False

INTERNET_ON = _internet()


def assert_something_resembling_minimization_happened(p0, e0, traj, mol):
    """ Raises assertion error if results do not correspond to a successful optimization

    Args:
        p0 (Array[length]): initial position array
        e0 (Scalar[energy]): initial energy
        traj (moldesign.Trajectory): trajectory created from minimization
        mol (moldesign.Molecule): state of molecule AFTER minimization

    Returns:

    """
    # Use assert_almost_equal to account for small numerical differences in
    # convergence noise, serialization methods, etc.

    import scipy.optimize.optimize

    assert traj.num_frames > 1

    assert_almost_equal(traj.potential_energy[0],
                        e0, decimal=9)
    assert_almost_equal(traj.potential_energy[-1],
                        mol.potential_energy, decimal=9)
    assert traj.potential_energy[-1] < e0

    assert (traj.positions[0] == p0).all()
    assert (traj.positions[-1] != p0).any()
    assert (traj.positions[-1] == mol.positions).all()

    recalc = mol.calculate_potential_energy(usecache=False)
    assert (recalc - traj.potential_energy[-1]) < 1.0e-8 * u.hartree

    scipyresult = getattr(traj, 'info', None)
    if isinstance(scipyresult, scipy.optimize.optimize.OptimizeResult):
        np.testing.assert_allclose(scipyresult.x,
                                   mol.positions.defunits_value().flat)
        assert_almost_equal(scipyresult.fun,
                            mol.potential_energy.defunits_value())

    if mol.constraints or mol.energy_model.params.get('constrain_hbonds', False):
        assert_constraints_satisfied(traj, p0, mol)


def assert_something_resembling_dynamics_happened(traj, mol, p0, t0, duration):
    """ Checks that dynamics routines have fulfilled their obligations.
    """
    assert mol is traj.mol
    assert traj.time[0] == t0
    assert mol.time == traj.time[-1]

    # Energy test is not exact due to some non-deterministic kernels (e.g., OpenMM CPU model)
    assert abs(traj.potential_energy[-1] - mol.calculate_potential_energy()) < 1e-10 * u.eV
    assert_almost_equal(traj.positions[-1], mol.positions)

    # lower precision for first step due to noise from single-precision positions (i.e. OpenMM CPU)
    assert_almost_equal(traj.positions[0], p0, decimal=5)

    lasttime = -np.inf * u.fs
    for step in traj:
        assert step.time > lasttime

    # currently we permit the integrator to be a little off in how long it integrates for
    assert mol.time - traj.time[0] > (duration - mol.integrator.params.frame_interval)

    if mol.constraints or mol.energy_model.params.get('constrain_hbonds', False):
        assert_constraints_satisfied(traj, p0, mol)


def assert_constraints_satisfied(traj, p0, mol):
    """ Checks that constraints were satisfied during molecular motion
    """
    for constraint in mol.constraints:
        assert constraint.satisfied()

    # TODO: check whole trajectory (not necessarily conserved during minimization however)
    # TODO: check hbonds, if energy model calls for it


def assert_almost_equal(actual, desired, **kwargs):
    units = mdt.units.get_units(actual)

    np.testing.assert_almost_equal(units.value_of(actual),
                                   units.value_of(desired),
                                   **kwargs)

