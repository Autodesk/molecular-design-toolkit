from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library
standard_library.install_aliases()

# Copyright 2017 Autodesk Inc.
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
import itertools

import numpy as np

import moldesign as mdt
from .. import units as u
from ..method import Method


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

    _CALLS_MDT_IN_DOCKER = False  # gets set to true if a python-interfaced dependency is missing

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
        If the energy model provides its own minimizer, it should be hooked up here
        """
        raise NotImplementedError()

    def get_formal_charge(self):
        """Determine the formal charge of the molecular system.
         This can be set either as a molecular attribute OR in the parameters of the energy model.

         Returns:
             u.Scalar[charge]: the formal charge used for this model
        """
        if 'charge' in self.params and self.params.charge is not None:
            return self.params.charge.value_in(u.q_e)
        elif hasattr(self.mol, 'charge'):
            return self.mol.charge.value_in(u.q_e)
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

        for iatom, idim in itertools.product(range(self.mol.num_atoms), range(3)):
            print('\rFinite differencing %s for atom %d/%d'%('xyz'[iatom],
                                                             iatom+1,
                                                             self.mol.num_atoms), end=' ')
            if direction == 0:
                self.mol.positions[iatom, idim] += stepsize / 2.0
                eplus = self.mol.calc_potential_energy()
                pplus = self.mol.properties

                self.mol.positions[iatom, idim] -= stepsize
                eminus = self.mol.calc_potential_energy()
                pminus = self.mol.properties

                self.mol.positions[iatom, idim] += stepsize / 2.0  # resets the position
                properties.append((pminus, pplus))
                forces[iatom, idim] = (eminus - eplus) / stepsize

            elif direction in (-1, 1):
                self.mol.positions[iatom, idim] += direction * stepsize
                enew = self.mol.calc_potential_energy()
                properties.append(self.mol.properties)
                forces[iatom, idim] = (e0 - enew) / (direction * stepsize)
                self.mol.positions[iatom, idim] -= direction * stepsize

        return forces, properties

    def prep(self):
        """
        Prepare to run. Possibly do a test to ensure that the model is ready.
        """
        raise NotImplementedError('EnergyModelBase is an abstract base class')


class MMBase(EnergyModelBase):
    """Common interface for molecular mechanics"""

    PARAMETERS = (EnergyModelBase.PARAMETERS +
                  list(mdt.parameters.mm_model_parameters.values()))

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mdtforcefield = None


class QMBase(EnergyModelBase):
    """Common interface for quantum mechanics"""

    PARAMETERS = list(mdt.parameters.qm_model_parameters.values())

    DEFAULT_PROPERTIES = ['potential_energy',
                          'nuclear_repulsion',
                          'dipole_moment',
                          'orbitals',
                          'orbital_energies']
    ALL_PROPERTIES = DEFAULT_PROPERTIES
    # properties will be a pretty long list for most packages

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


