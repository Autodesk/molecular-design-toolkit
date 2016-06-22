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
import scipy.optimize

import moldesign as mdt
from moldesign import units as u


class MinimizerBase(object):

    _strip_units = True  # callbacks expect and return dimensionless quantities scaled to default unit system
    constraint_restraints = True  # if True, add restraint penalties for both constraints and restraints

    def __init__(self, mol, nsteps=20,
                 force_tolerance=moldesign.core.data.DEFAULT_FORCE_TOLERANCE,
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


class BFGS(MinimizerBase):
    _strip_units = True

    def run(self):
        print 'Starting geometry optimization: scipy.fmin_bfgs with %s gradients' % self.gradtype
        options = {'disp': True}
        if self.nsteps is not None:
            options['maxiter'] = self.nsteps

        if self.gradtype == 'analytical': grad = self.grad
        else: grad = None

        if self.force_tolerance is not None:
            options['gtol'] = self.force_tolerance.defunits().magnitude

        result = scipy.optimize.minimize(self.objective,
                                         self.mol.positions.defunits().magnitude,
                                         method='bfgs',
                                         jac=grad,
                                         callback=self.callback,
                                         options=options)

        options['maxiter'] = max(options['maxiter']/5, 2)  # reduce number of steps if we need to enforce constraints
        for i in xrange(3):  # try to satisfy constraints
            if all(c.satisfied() for c in self.mol.constraints): break
            self.restraint_multiplier += 0.5
            print 'Constraints not satisfied - new restraint multiplier %f ...' % \
                  self.restraint_multiplier
            result = scipy.optimize.minimize(self.objective,
                                 self.mol.positions.defunits().magnitude,
                                 method='bfgs',
                                 jac=grad,
                                 callback=self.callback,
                                 options=options)
        if not all(c.satisfied() for c in self.mol.constraints):
            print 'WARNING! Constraints not satisfied'
        self.traj.info = result


class GradientDescent(MinimizerBase):
    _strip_units = False

    def __init__(self, mol, max_atom_move=0.075 * u.angstrom,
                 gamma=0.05 * (u.angstrom ** 2) / u.eV,
                 **kwargs):
        """
        :param max_atom_move: Maximum amount to move an atom
        """
        super(GradientDescent, self).__init__(mol, **kwargs)
        assert 'forces' in self.request_list, 'Gradient descent built-in gradients'
        self.max_atom_move = max_atom_move
        self.gamma = gamma
        self._last_energy = None

    def run(self):
        self.run_once()
        self.nsteps = max(self.nsteps / 5, 2)
        for i in xrange(3):  # try to satisfy constraints
            if all(c.satisfied() for c in self.mol.constraints): break
            self.restraint_multiplier += 0.5
            print '\nWarning: constraints not satisfied - new restraint multiplier %f ...' % \
                  self.restraint_multiplier
            self.run_once()

        if not all(c.satisfied() for c in self.mol.constraints):
            print 'WARNING! Constraints not satisfied'

    def run_once(self):
        print 'Starting geometry optimization: built-in gradient descent'
        laststep = self.mol.calculate()
        lastenergy = self.objective(self.mol.positions)
        for i in xrange(self.nsteps):
            grad = self.grad(self.mol.positions)
            if np.abs(grad.max()) < self.force_tolerance:
                return
            move = -self.gamma * grad
            mmax = np.abs(move).max()
            if mmax > self.max_atom_move:  # rescale the move
                scale = self.max_atom_move / mmax
                print 'Move too big: scaling by step by {scale.magnitude:.6f}'.format(scale=scale)
                move *= scale
            self.mol.positions += move
            newenergy = self.objective(self.mol.positions)
            if newenergy > lastenergy:
                print 'Energy increased by {x.magnitude:.3e} {x.units}!'.format(x=(newenergy-lastenergy)) + \
                      ' Reverting to previous step and reducing step size.'
                self.gamma /= 2.0
                self.max_atom_move /= 2.0
                self.mol.positions = laststep['positions']
                self.mol.properties = laststep
            else:
                laststep = self.mol.calculate(self.request_list)
                lastenergy = newenergy
            self.callback()



def bfgs(*args, **kwargs):
    """
    The function version of BFGS
    """
    minobj = BFGS(*args, **kwargs)
    return minobj()


def gradient_descent(*args, **kwargs):
    minobj = GradientDescent(*args, **kwargs)
    return minobj()


lookup = {'gradient descent': gradient_descent,
          'bfgs': bfgs}
