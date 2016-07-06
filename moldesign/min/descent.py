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

from moldesign import utils
from moldesign import units as u
from .base import MinimizerBase
from . import toplevel


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class GradientDescent(MinimizerBase):
    """ A careful (perhaps overly careful) gradient descent implementation designed to relax
    structures far from equilibrium.

    The structure is adjusted along the force direction, but the step size is capped by
        ``max_atom_move`` (.075 angstrom by default). If the energy increases after a move,
        the move is repeated with half the step size.

    Note:
        This algorithm is good at stably removing large forces, but it's very poorly suited to
           locating any type of critical point.
    """

    _strip_units = False

    @utils.args_from(MinimizerBase,
                       inject_kwargs={'max_atom_move': 0.075*u.angstrom,
                                      'gamma': 0.05*(u.angstrom**2)/u.eV})
    def __init__(self, mol, max_atom_move=0.075*u.angstrom,
                 gamma=0.05*(u.angstrom ** 2)/u.eV,
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
        lastproperties = self.mol.calculate()
        lastenergy = self.objective(self._coords_to_vector(self.mol.positions))
        current = self._coords_to_vector(self.mol.positions)

        for i in xrange(self.nsteps):
            grad = self.grad(current)
            if np.abs(grad.max()) < self.force_tolerance:
                return
            move = -self.gamma * grad
            mmax = np.abs(move).max()
            if mmax > self.max_atom_move:  # rescale the move
                scale = self.max_atom_move / mmax
                #print 'Move too big: scaling step by {scale.magnitude:.6f}'.format(scale=scale)
                move *= scale
            current += move
            newenergy = self.objective(current)
            if newenergy > lastenergy:
                print 'Energy increased by {x.magnitude:.3e} {x.units}!'.format(
                        x=(newenergy-lastenergy)) + \
                      ' Reverting to previous step and reducing step size.'
                self.gamma /= 2.0
                self.max_atom_move /= 2.0
                self.mol.positions = lastproperties.positions
                self.mol.properties = lastproperties
            else:  # update with new positions
                self._sync_positions(current)
                lastproperties = self.mol.calculate(self.request_list)
                lastenergy = newenergy
            self.callback()


gradient_descent = GradientDescent._as_function('gradient_descent')
exports(gradient_descent)
toplevel(gradient_descent)

