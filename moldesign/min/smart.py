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
        super().__init__(*args, **kwargs)

    def run(self):
        # If forces are already low, go directly to the quadratic convergence methods and return
        forces = self.mol.calculate_forces()
        if abs(forces).max() <= self.gd_threshold:
            spmin = self._make_quadratic_method()
            spmin.run()
            self.traj = spmin.traj
            self.current_step = spmin.current_step
            return

        # Otherwise, remove large forces with gradient descent; exit if we pass the cycle limit
        descent_kwargs = self.kwargs.copy()
        descent_kwargs['force_tolerance'] = self.gd_threshold
        descender = GradientDescent(*self.args, **descent_kwargs)
        descender.run()
        if descender.current_step >= self.nsteps:
            self.traj = descender.traj
            return

        # Finally, use a quadratic method to converge the optimization
        kwargs = dict(_restart_from=descender.current_step,
                      _restart_energy=descender._initial_energy)
        kwargs['frame_interval'] = self.kwargs.get('frame_interval',
                                                   descender.frame_interval)
        spmin = self._make_quadratic_method(kwargs)
        spmin.current_step = descender.current_step
        spmin.run()
        self.traj = descender.traj + spmin.traj
        self.current_step = spmin.current_step

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

