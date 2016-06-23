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
import sys

import numpy as np

import moldesign as mdt
from moldesign import data
from moldesign import units as u


class MinimizerBase(object):

    _strip_units = True  # callbacks expect and return dimensionless quantities scaled to default unit system
    constraint_restraints = True  # if True, add restraint penalties for both constraints and restraints

    def __init__(self, mol, nsteps=20,
                 force_tolerance = data.DEFAULT_FORCE_TOLERANCE,
                 frame_interval=None,
                 restraint_multiplier=1.0):
        self.mol = mol
        self.nsteps = nsteps
        self.force_tolerance = force_tolerance
        self.frame_interval = mdt.utils.if_not_none(frame_interval, nsteps/10)

        # Set up the trajectory to track the minimization
        self.traj = mdt.Trajectory(mol)
        self.current_step = 0
        self.traj.new_frame(minimization_step=0,
                            annotation='minimization steps:%d (energy=%s)' %
                                       (0, mol.calc_potential_energy()))
        self._initial_energy = None
        self._last_energy = None
        self._last_grad = None

        # Figure out whether we'll use gradients
        self.request_list = ['potential_energy']
        self.restraint_multiplier = restraint_multiplier
        if 'forces' in mol.energy_model.ALL_PROPERTIES:
            self.gradtype = 'analytical'
            self.request_list.append('forces')
        else:
            self.gradtype = 'approximate'
            assert len(mol.constraints) == 0, \
                'Constrained minimization only available with analytical gradients'

    def objective(self, coords):
        if self._strip_units:
            self.mol.positions = coords * u.default.length
        else:
            self.mol.positions = coords
        self.mol.calculate(requests=self.request_list)
        pot = self.mol.potential_energy

        if self.constraint_restraints:
            for constraint in self.mol.constraints:
                pot += self.restraint_multiplier * constraint.restraint_penalty()

        if self._initial_energy is None: self._initial_energy = pot
        self._last_energy = pot
        if self._strip_units: return pot.defunits().magnitude
        else: return pot.defunits()

    def grad(self, coords):
        if self._strip_units:
            self.mol.positions = coords * u.default.length
        else:
            self.mol.positions = coords
        self.mol.calculate(requests=self.request_list)
        grad = -self.mol.forces

        if self.constraint_restraints:
            for constraint in self.mol.constraints:
                grad -= self.restraint_multiplier * constraint.restraint_penalty_force()

        self._last_grad = grad
        if self._strip_units: return grad.defunits().magnitude
        else: return grad.defunits()

    def __call__(self):
        self.run()
        self.traj.new_frame(minimization_step=self.current_step,
               annotation='minimization result (%d steps) (energy=%s)' %
                          (self.current_step, self.mol.potential_energy))
        return self.traj

    def run(self):
        raise NotImplementedError('This is an abstract base class')

    def callback(self, *args):
        self.current_step += 1
        if self.current_step % self.frame_interval != 0: return

        self.mol.calculate(self.request_list)
        self.traj.new_frame(minimization_step=self.current_step,
                       annotation='minimization steps:%d (energy=%s)' %
                                  (self.current_step, self.mol.potential_energy))
        if self.nsteps is None:
            message = ['Minimization step %d' % self.current_step]
        else:
            message = ['Step %d/%d' % (self.current_step, self.nsteps)]

        if self._last_energy is not None:
            message.append(u'\u0394E={x.magnitude:.3e} {x.units}'.format(
                x=self._last_energy - self._initial_energy))

        if self.gradtype == 'analytical' and self._last_grad is not None:
            force = self._last_grad

            message.append(u'RMS \u2207E={rmsf.magnitude:.3e}, '
                           u'max \u2207E={mf.magnitude:.3e}{mf.units}'.format(
                rmsf=np.sqrt(force.dot(force) / self.mol.ndims),
                mf=np.abs(force).max()))

        if self.constraint_restraints and self.mol.constraints:
            nsatisfied = 0
            for c in self.mol.constraints:
                if c.satisfied(): nsatisfied += 1
            message.append('constraints:%d/%d' % (nsatisfied, len(self.mol.constraints)))

        print ', '.join(message)
        sys.stdout.flush()
