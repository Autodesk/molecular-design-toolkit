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

from .base import MinimizerBase
from . import toplevel, BFGS, GradientDescent, SequentialLeastSquares

from moldesign import utils
from moldesign import units as u
from moldesign.utils import exports

GDTHRESH = 1.5*u.eV/u.angstrom

@exports
class SmartMin(MinimizerBase):
    """ Uses gradient descent until forces fall below a threshold,
    then switches to BFGS (unconstrained) or SLSQP (constrained).

    Args:
        gd_threshold (u.Scalar[force]): Use gradient descent if there are any forces larger
           than this; use an approximate hessian method (BFGS or SLSQP) otherwise

    Note:
        Not really that smart.
    """
    # TODO: use non-gradient methods if forces aren't available

    _strip_units = True

    @utils.args_from(MinimizerBase,
                     inject_kwargs={'gd_threshold': GDTHRESH})
    def __init__(self, *args, **kwargs):
        self.gd_threshold = kwargs.pop('gd_threshold', GDTHRESH)
        self.args = args
        self.kwargs = kwargs
        self._spmin = None
        self._descent = None
        self._traj = None
        self._currentstep = None
        self.__foundmin = None
        super().__init__(*args, **kwargs)

    def _run(self):
        # If forces are already low, go directly to the quadratic convergence methods and return
        forces = self.mol.calculate_forces()
        if abs(forces).max() <= self.gd_threshold:
            self._spmin = self._make_quadratic_method()
            self._spmin._run()
            self.traj = self._spmin.traj
            return

        # Otherwise, remove large forces with gradient descent; exit if we pass the cycle limit
        descent_kwargs = self.kwargs.copy()
        descent_kwargs['force_tolerance'] = self.gd_threshold
        self._descender = GradientDescent(*self.args, **descent_kwargs)
        self._descender._run()
        if self._descender.current_step >= self.nsteps:
            self.traj = self._descender.traj
            return

        # Finally, use a quadratic method to converge the optimization
        kwargs = dict(_restart_from=self._descender.current_step,
                      _restart_energy=self._descender._initial_energy)
        kwargs['frame_interval'] = self.kwargs.get('frame_interval',
                                                   self._descender.frame_interval)
        self._spmin = self._make_quadratic_method(kwargs)
        self._spmin.current_step = self.current_step
        self._spmin._foundmin = self._foundmin
        self._spmin._run()
        self.traj = self._descender.traj + self._spmin.traj
        self.traj.info = getattr(self._spmin, 'info', None)

    @property
    def _foundmin(self):
        if self._descent:
            return self._descent._foundmin
        elif self._spmin:
            return self._spmin._foundmin
        else:
            return self.__foundmin

    @_foundmin.setter
    def _foundmin(self, val):
        self.__foundmin = val

    @property
    def currentstep(self):
        if self._descent:
            return self._descent.currentstep
        elif self._spmin:
            return self._spmin.currentstep
        else:
            return self._currentstep

    @currentstep.setter
    def currentstep(self, val):
        self._currentstep = val

    def _make_quadratic_method(self, kwargs=None):
        if kwargs is None: kwargs = {}
        kw = self.kwargs.copy()
        kw.update(kwargs)
        if self.mol.constraints:
            spmin = SequentialLeastSquares(*self.args, **kw)
        else:
            spmin = BFGS(*self.args, **kw)
        return spmin


minimize = SmartMin._as_function('minimize')
exports(minimize)
toplevel(minimize)

