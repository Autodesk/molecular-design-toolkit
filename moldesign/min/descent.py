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
import numpy as np

import moldesign as mdt
from .. import units as u
from ..utils import exports
from .base import MinimizerBase
from . import toplevel


@exports
class GradientDescent(MinimizerBase):
    """ A careful (perhaps overly careful) gradient descent implementation designed to relax
    structures far from equilibrium.

    A backtracking line search is performed along the steepest gradient direction.

    The maximum move for any single atom is also limited by ``max_atom_move``

    Note:
        This algorithm is good at stably removing large forces, but it's very poorly suited to
           locating any type of critical point; don't use this to find a minimum!

    References:
        https://www.math.washington.edu/~burke/crs/408/lectures/L7-line-search.pdf

    Args:
        mol (moldesign.Molecule): molecule to minimize
        max_atom_move (Scalar[length]): maximum displacement of a single atom
        scaling (Scalar[length/force]): unit of displacement per unit force
        gamma (float): number between 0 and 1 indicating scale factor for backtracking search
        control (float): threshold for terminating line search; this is a proportion
             (0<=``control``<=1) of the expected function decrease
        **kwargs (dict): kwargs from :class:`MinimizerBase`
    """

    _strip_units = False

    def __init__(self, mol,
                 max_atom_move=0.05*u.angstrom,
                 scaling=0.01*u.angstrom**2/u.eV,
                 gamma=0.4, control=0.25,
                 **kwargs):
        super().__init__(mol, **kwargs)
        assert 'forces' in self.request_list, 'Gradient descent built-in gradients'
        self.max_atom_move = max_atom_move
        self.scaling = scaling
        self.gamma = gamma
        self.control = control
        self._last_energy = None

    def _run(self):
        print('Starting geometry optimization: built-in gradient descent')
        lastenergy = self.objective(self._coords_to_vector(self.mol.positions))
        current = self._coords_to_vector(self.mol.positions)

        for i in range(self.nsteps):
            grad = self.grad(current)
            if np.abs(grad.max()) < self.force_tolerance:  # converged
                return

            move = self.scale_move(grad)
            armijo_goldstein_prefac = self.control * move.norm()

            for icycle in range(0, 10):
                g = self.gamma**icycle
                newpos = self._make_move(current, g * move)

                # move direction may be different than gradient direction due to constraints
                move_vec = (newpos-current).normalized()
                if grad.dot(move_vec) >= 0.0:  # move flipped direction!
                    if self._constraint_convergence(newpos, current, grad):
                        return  # flip was because we're converged
                    else:  # flip was because move was too big
                        newenergy = np.inf * u.default.energy
                        continue

                try:
                    newenergy = self.objective(newpos)
                except mdt.QMConvergenceError:
                    continue

                if newenergy <= lastenergy + g * armijo_goldstein_prefac * grad.dot(move_vec):
                    break
            else:
                if newenergy >= lastenergy:
                    raise mdt.ConvergenceFailure('Line search failed')

            if self._constraint_convergence(newpos, current, grad):
                return
            else:
                current = newpos
                lastenergy = newenergy
                self._sync_positions(current)

            self.callback()

    def scale_move(self, grad):
        move = -self.scaling*grad
        mmax = np.abs(move).max()
        if mmax > self.max_atom_move:  # rescale the move
            move *= self.max_atom_move/mmax
        return move

    def _make_move(self, current, move):
        if self.mol.constraints:
            # TODO: get constraint forces from lagrange multipliers and use them to check for convergence
            self._sync_positions(current)
            prev = self.mol.positions.copy()
            self._sync_positions(current+move)

            mdt.geom.shake_positions(self.mol, prev)
            return self._coords_to_vector(self.mol.positions)
        else:
            return current + move

    def _constraint_convergence(self, pos, lastpos, energygrad):
        """ Test for force-based convergence after projecting out constraint forces

        Until the shake method starts explicitly storing constraint forces, we calculate this
        direction as the SHAKE-adjusted displacement vector from the current descent step
        """
        direction = mdt.mathutils.normalized((pos - lastpos).flatten())
        proj_grad = energygrad.dot(direction)
        return abs(proj_grad) < self.force_tolerance


gradient_descent = GradientDescent._as_function('gradient_descent')
exports(gradient_descent)
toplevel(gradient_descent)

