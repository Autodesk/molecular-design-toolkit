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

from .. import units as u
from ..utils import exports
from .base import MinimizerBase
from .. import exceptions
from . import toplevel


class ScipyMinimizer(MinimizerBase):
    """ SciPy's implementation of the BFGS method, with gradients if available.

    Args:
        bfgs_threshold (u.Scalar[force]): Maximum force on a single atom

    Note:
        This implementation will fail rapidly if large forces are present (>> 1 eV/angstrom).
    """
    _strip_units = True
    _METHOD_NAME = None
    _TAKES_FTOL = False
    _TAKES_GTOL = False

    def _run(self):
        import scipy.optimize

        if self.mol.constraints and self._METHOD_NAME == 'bfgs':
            raise exceptions.NotSupportedError('BFGS minimization does not '
                                               'support constrained minimization')

        print('Starting geometry optimization: SciPy/%s with %s gradients'%(
            self._METHOD_NAME, self.gradtype))
        options = {'disp': True}
        if self.nsteps is not None:
            options['maxiter'] = self.nsteps

        if self.gradtype == 'analytical':
            grad = self.grad
        else: grad = None

        if self.force_tolerance is not None:
            if self._TAKES_GTOL:
                options['gtol'] = self.force_tolerance.defunits().magnitude
            elif self._TAKES_FTOL:
                print('WARNING: this method does not use force to measure convergence; '
                      'approximating force_tolerance keyword')
                options['ftol'] = (self.force_tolerance * u.angstrom / 50.0).defunits_value()
            else:
                print('WARNING: no convergence criteria for this method; using defaults')

        self._optimize_kwargs = dict(method=self._METHOD_NAME,
                                     options=options)
        self._constraint_multiplier = 1.0

        result = scipy.optimize.minimize(self.objective,
                                         self._coords_to_vector(self.mol.positions),
                                         jac=grad,
                                         callback=self.callback,
                                         constraints=self._make_constraints(),
                                         **self._optimize_kwargs)

        if self.mol.constraints:
            result = self._force_constraint_convergence(result)

        self.traj.info = result

        finalprops = self._calc_cache[tuple(result.x)]
        self.mol.positions = finalprops.positions
        self.mol.properties = finalprops

    def _force_constraint_convergence(self, result):
        """ Make sure that all constraints are satisfied, ramp up the constraint functions if not

        Note - if additional iterations are necessary, this will destroy the scipy optimize results
        object stored at self.traj.info. Not sure what to do about that
        """
        import scipy.optimize

        for i in range(5):
            for constraint in self.mol.constraints:
                if not constraint.satisfied():
                    break
            else:
                return result

            print('Constraints not satisfied; raising penalties ...')

            self._constraint_multiplier *= 10.0
            result = scipy.optimize.minimize(self.objective,
                                             self._coords_to_vector(self.mol.positions),
                                             jac=self.grad if self.gradtype=='analytical' else None,
                                             callback=self.callback,
                                             constraints=self._make_constraints(),
                                             **self._optimize_kwargs)
        return result

    def _make_constraints(self):
        from .. import geom

        constraints = []
        for constraint in geom.get_base_constraints(self.mol.constraints):
            fun, jac = self._make_constraint_funs(constraint)
            constraints.append(dict(type='eq',
                                    fun=fun,
                                    jac=jac))
        return constraints


    def _make_constraint_funs(self, const):
        def fun(v):
            self._sync_positions(v)
            return const.error().defunits_value() * self._constraint_multiplier

        def jac(v):
            self._sync_positions(v)
            return (const.gradient().defunits_value().reshape(self.mol.num_atoms*3)
                    * self._constraint_multiplier)

        return fun, jac


@exports
class BFGS(ScipyMinimizer):
    _METHOD_NAME = 'bfgs'
    _TAKES_GTOL = True


bfgs = BFGS._as_function('bfgs')
exports(bfgs)
toplevel(bfgs)


@exports
class SequentialLeastSquares(ScipyMinimizer):
    _METHOD_NAME = 'SLSQP'
    _TAKES_FTOL = True


sequential_least_squares = SequentialLeastSquares._as_function('sequential_least_squares')
exports(sequential_least_squares)
toplevel(sequential_least_squares)
