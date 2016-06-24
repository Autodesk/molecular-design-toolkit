# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import numpy as np

import moldesign as mdt
from moldesign import units as u
from moldesign.keywords import mm_model_parameters as mmp, qm_model_parameters as qmp
from moldesign.method import Method


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

    PARAMETERS = EnergyModelBase.PARAMETERS+[
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


