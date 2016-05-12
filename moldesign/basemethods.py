"""
This module contains abstract base classes for potential models, integrators, and various
associated data types (force fields, orbitals, basis sets, etc.).
"""
import sys

import numpy as np
import moldesign as mdt
from moldesign import units as u
from moldesign.keywords import mm_model_parameters as mmp
from moldesign.keywords import qm_model_parameters as qmp
from moldesign.keywords import integrator_parameters as ips
from moldesign.utils import DotDict, if_not_none


def _make_method(cls, params, mol):
    """
    Helper for __reduce__ to use kwargs
    """
    obj = cls(**params)
    obj.mol = mol
    return obj


class Method(object):
    """Base class for energy models and integrators"""
    PARAMETERS = DotDict()
    """
    A list of Parameter objects describing the keyword arguments that this class
    recognizes.

    Attributes:
       mol (mdt.Molecule): the molecule this method is associated with
    """

    def __reduce__(self):
        return _make_method, (self.__class__, self.params, self.mol)

    def __init__(self, **params):
        """
        :param params:
        :return:
        """
        # TODO: better documentation for the expected keywords

        self._prepped = False
        self.status = None
        self.mol = None
        self.params = DotDict(params)
        # Set default parameter values
        for param in self.PARAMETERS:
            if param.name not in self.params:
                self.params[param.name] = param.default

    @classmethod
    def get_parameters(cls):
        """
        This doesn't do anything right now except provide guidelines for programmers
        """
        return cls.PARAMETERS

    def get_forcefield(self):
        raise NotImplementedError()

    @classmethod
    def print_parameters(cls):
        params = cls.PARAMETERS
        lines = []
        for obj in params:
            description = ''
            if obj.choices:
                description = '%s' % obj.choices
                if obj.types:
                    description += ' or '
            if obj.types:
                description += 'Type %s' % obj.types

            doc = '%s: %s (DEFAULT: %s)' % (obj.name, description, obj.default)
            lines.append(doc)
        return '\n'.join(lines)


####################################
# The first kind of Method is an "EnergyModel" - something that calculates
# potential energy (and other properties, usually) as a function of atomic positions.
class EnergyModelBase(Method):
    """
    Base class for all energy models
    """

    # TODO:should add some architecture to check implementations, e.g.
    # TODO: ensure prep gets called when necessary, make sure all parameters can be consumed

    DEFAULT_PROPERTIES = ['potential_energy', 'forces']
    """List[str]: list of the properties that are always calculated by this method"""
    ALL_PROPERTIES = DEFAULT_PROPERTIES
    """List[str]: List of all the properties that this model can calculate"""

    PARAMETERS = []

    def calculate(self, requests):
        """Calculate the the default properties and any additiona requests

        Arguments:
            requests (List[str]): the requested properties to calculate

        Returns:
           utils.DotDict: A dict of calculated properties (or a job object that will return them)
        """
        self.prep()
        raise NotImplementedError('EnergyModelBase is an abstract base class')

    def minimize(self, method='bfgs', **kwargs):
        """
        Relax the parent molecule's energy a built-in minimization scheme.
        Many energy methods will provide their own minimization scheme and will override this method
        """
        if method == 'bfgs':
            return mdt.minimizers.bfgs(self.mol, **kwargs)
        elif method == 'gradient descent':
            return mdt.minimizers.gradient_descent(self.mol, **kwargs)
        else:
            raise ValueError('Unknown minimization method %s' % method)

    def get_formal_charge(self):
        """Determine the formal charge of the molecular system.
         This can be set either as a molecular attribute OR in the parameters of the energy model.

         Returns:
             u.Scalar[charge]: the formal charge used for this model
        """
        if 'charge' in self.params and self.params.charge is not None:
            return self.params.charge
        elif hasattr(self.mol, 'charge'):
            return self.mol.charge
        else:
            return 0

    def finite_difference_force(self, direction=0, stepsize=0.025 * u.angstrom):
        """
        Compute force using a finite difference with the given step size.

        Args:
            direction (int): EITHER +1, -1, (for one-sided finite differences) or
               0 (for central  difference - better but twice as expensive)
            step (u.Scalar[lenght]): step size to take in each direction

        Returns:
            u.Vector[force]: force vector, len= `self.mol.ndims`
        """
        # TODO: this should totally be parallelized - how do we request/configure/control this?
        forces = np.zeros(self.mol.ndims) * u.default.force
        properties = []
        if direction in (-1, 1):
            e0 = self.mol.calc_potential_energy()
        else:
            assert direction == 0, 'Finite difference direction must be -1, 0, or 1'

        for idim in xrange(self.mol.ndims):
            print '\rFinite difference step %d/%d' % (idim + 1, self.mol.ndims),
            if direction == 0:
                self.mol.positions[idim] += stepsize / 2.0
                eplus = self.mol.calc_potential_energy()
                pplus = self.mol.properties

                self.mol.positions[idim] -= stepsize
                eminus = self.mol.calc_potential_energy()
                pminus = self.mol.properties

                self.mol.positions[idim] += stepsize / 2.0  # resets the position
                properties.append((pminus, pplus))
                forces[idim] = (eminus - eplus) / stepsize

            elif direction in (-1, 1):
                self.mol.positions[idim] += direction * stepsize
                enew = self.mol.calc_potential_energy()
                properties.append(self.mol.properties)
                forces[idim] = (e0 - enew) / (direction * stepsize)
                self.mol.positions[idim] -= direction * stepsize

        return forces, properties

    def prep(self):
        """
        Prepare to run. Possibly do a test to ensure that the model is ready.
        """
        raise NotImplementedError('EnergyModelBase is an abstract base class')


class MMBase(EnergyModelBase):
    """Common interface for molecular mechanics"""

    PARAMETERS = EnergyModelBase.PARAMETERS + [
        mmp.forcefield, mmp.implicit_solvent,
        mmp.cutoff, mmp.nonbonded, mmp.constrain_hbonds,
        mmp.constrain_water,
        mmp.solute_dielectric, mmp.solvent_dielectric,
        mmp.ewald_error, mmp.periodic]

    def __init__(self, *args, **kwargs):
        super(MMBase, self).__init__(*args, **kwargs)
        self.mdtforcefield = None


class QMBase(EnergyModelBase):
    """Common interface for quantum mechanics"""

    DEFAULT_PROPERTIES = ['potential_energy',
                          'nuclear_repulsion',
                          'dipole_moment',
                          'orbitals',
                          'orbital_energies']
    ALL_PROPERTIES = DEFAULT_PROPERTIES
    # properties will be a pretty long list for most packages

    PARAMETERS = [qmp.theory, qmp.charge, qmp.multiplicity,
                  qmp.basis_set, qmp.symmetry, qmp.wfn_guess]

    def set_wfn_guess(self):
        raise NotImplementedError


class QMMMBase(EnergyModelBase):
    DEFAULT_PROPERTIES = ['potential_energy',
                          'qm_energy',
                          'mm_energy',
                          'interaction_energy'
                          'qm_dipole_moment',
                          'orbitals',
                          'orbital_energies']
    ALL_PROPERTIES = DEFAULT_PROPERTIES

    PARAMETERS = MMBase.PARAMETERS + QMBase.PARAMETERS



#################################
# Next are some base integrators
class IntegratorBase(Method):
    """Base class for all integrators"""

    PARAMETERS = [ips.timestep, ips.frame_interval]

    def run(self, run_for):
        """
        To be called by parent molecule
        :param run_for: number of steps (if integer), or amount of time (if has units of time)
        :return: trajectory
        """
        raise NotImplementedError('This is an abstract base class!')

    def prep(self):
        """
        Prepare to run. Possibly do a test run.
        This might need to call the mol.model.build. Make sure you don't have a
        circular call here
        :return:
        """
        raise NotImplementedError()

    @staticmethod
    def time_to_steps(time, timestep):
        try:
            dims = time.dimensionality
            assert len(dims) == 1 and dims['[time]'] == 1.0
        except (AttributeError, AssertionError):
            assert type(time) == int, "argument to integrator.run must have units of time or be an int"
            return time
        else:
            return int(time / timestep)


class MDBase(IntegratorBase):
    PARAMETERS = (IntegratorBase.PARAMETERS +
                  [ips.remove_translation,ips.remove_rotation])


class ConstantTemperatureBase(MDBase):
    PARAMETERS = MDBase.PARAMETERS + [ips.temperature]


class LangevinBase(ConstantTemperatureBase):
    PARAMETERS = ConstantTemperatureBase.PARAMETERS + [ips.collision_rate]

