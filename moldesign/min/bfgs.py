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
import scipy.optimize

from .base import MinimizerBase
from . import toplevel
from moldesign import utils


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class BFGS(MinimizerBase):
    """ SciPy's implementation of the BFGS method, with gradients if available.

    Args:
        bfgs_threshold (u.Scalar[force]): Maximum force on a single atom

    Note:
        This implementation will fail rapidly if large forces are present (>> 1 eV/angstrom).
    """
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
                                         self._coords_to_vector(self.mol.positions),
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
                                 self._coords_to_vector(self.mol.positions),
                                 method='bfgs',
                                 jac=grad,
                                 callback=self.callback,
                                 options=options)
        if not all(c.satisfied() for c in self.mol.constraints):
            print 'WARNING! Constraints not satisfied'
        self.traj.info = result


bfgs = BFGS._as_function('bfgs')
exports(bfgs)
toplevel(bfgs)
