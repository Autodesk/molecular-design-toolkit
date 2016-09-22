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

    def __init__(self, mol, nsteps=20,
                 force_tolerance=data.DEFAULT_FORCE_TOLERANCE,
                 frame_interval=None,
                 _restart_from=0,
                 _restart_energy=None):
        self.mol = mol
        self.nsteps = nsteps - _restart_from
        self.force_tolerance = force_tolerance
        self.frame_interval = mdt.utils.if_not_none(frame_interval,
                                                    max(nsteps/10, 1))
        self._restart_from = _restart_from

        # Set up the trajectory to track the minimization
        self.traj = mdt.Trajectory(mol)
        self.current_step = _restart_from
        if _restart_energy is None:
            self.traj.new_frame(minimization_step=0,
                                annotation='minimization steps:%d (energy=%s)' %
                                           (0, mol.calc_potential_energy()))
        self._initial_energy = _restart_energy
        self._last_energy = None
        self._last_grad = None

        # Figure out whether we'll use gradients
        self.request_list = ['potential_energy']
        if 'forces' in mol.energy_model.ALL_PROPERTIES:
            self.gradtype = 'analytical'
            self.request_list.append('forces')
        else:
            self.gradtype = 'approximate'
            assert len(mol.constraints) == 0, \
                'Constrained minimization only available with analytical gradients'

    def _sync_positions(self, vector):
        """ Set the molecule's position
        """
        c = vector.reshape((self.mol.num_atoms, 3))
        if self._strip_units:
            self.mol.positions = c*u.default.length
        else:
            self.mol.positions = c

    def _coords_to_vector(self, coords):
        """ Convert position array to a flat vector
        """
        vec = coords.reshape(self.mol.num_atoms * 3).copy()
        if self._strip_units:
            return vec.magnitude
        else:
            return vec

    def objective(self, vector):
        """ Callback function to calculate the objective function
        """
        self._sync_positions(vector)
        try:
            self.mol.calculate(requests=self.request_list)
        except mdt.QMConvergenceError:  # returning infinity can help rescue some line searches
            return np.inf

        pot = self.mol.potential_energy

        if self._initial_energy is None: self._initial_energy = pot
        self._last_energy = pot
        if self._strip_units: return pot.defunits().magnitude
        else: return pot.defunits()

    def grad(self, vector):
        """ Callback function to calculate the objective's gradient
        """
        self._sync_positions(vector)
        self.mol.calculate(requests=self.request_list)
        grad = -self.mol.forces

        grad = grad.reshape(self.mol.num_atoms * 3)
        self._last_grad = grad
        if self._strip_units:
            return grad.defunits().magnitude
        else:
            return grad.defunits()

    def __call__(self):
        """ Run the minimization

        Returns:
            moldesign.Trajectory: the minimization trajectory
        """
        self.run()
        if self.traj.num_frames == 0 or self.traj[-1].minimization_step != self.current_step:
            self.traj.new_frame(minimization_step=self.current_step,
                                annotation='minimization result (%d steps) (energy=%s)'%
                                           (self.current_step, self.mol.potential_energy))
        return self.traj

    def run(self):
        raise NotImplementedError('This is an abstract base class')

    def callback(self, *args):
        """ To be called after each minimization step

        Args:
            *args: ignored
        """
        self.current_step += 1
        if self.current_step % self.frame_interval != 0: return

        self.mol.calculate(self.request_list)
        self.traj.new_frame(minimization_step=self.current_step,
                            annotation='minimization steps:%d (energy=%s)'%
                                       (self.current_step, self.mol.potential_energy))
        if self.nsteps is None:
            message = ['Minimization step %d' % self.current_step]
        else:
            message = ['Step %d/%d' % (self.current_step, self.nsteps + self._restart_from)]

        if self._last_energy is not None:
            message.append(u'\u0394E={x.magnitude:.3e} {x.units}'.format(
                x=self._last_energy - self._initial_energy))

        if self.gradtype == 'analytical' and self._last_grad is not None:
            force = self._last_grad

            message.append(u'RMS \u2207E={rmsf.magnitude:.3e}, '
                           u'max \u2207E={mf.magnitude:.3e} {mf.units}'.format(
                rmsf=np.sqrt(force.dot(force) / self.mol.ndims),
                mf=np.abs(force).max()))

        if self.mol.constraints:
            nsatisfied = 0
            for c in self.mol.constraints:
                if c.satisfied(): nsatisfied += 1
            message.append('constraints:%d/%d' % (nsatisfied, len(self.mol.constraints)))

        print(', '.join(message))
        sys.stdout.flush()

    @classmethod
    def _as_function(cls, newname):
        """ Create a function that runs this minimization
        """
        @mdt.utils.args_from(cls, allexcept=['self'])
        def asfn(*args, **kwargs):
            obj = cls(*args, **kwargs)
            return obj()
        asfn.__name__ = newname
        asfn.__doc__ = cls.__doc__

        return asfn




