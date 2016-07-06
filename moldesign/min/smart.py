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

from .base import MinimizerBase
from . import toplevel, BFGS, GradientDescent

from moldesign import utils
from moldesign import units as u

BFGSTHRESH = 1.5 * u.eV / u.angstrom

def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class SmartMin(MinimizerBase):
    """ Uses gradient descent until forces fall below a threshold,
    then switches to BFGS.

    Args:
        bfgs_threshold (u.Scalar[force]): Maximum force on a single atom

    Note:
        Not really that smart.
    """

    _strip_units = True

    @utils.args_from(MinimizerBase,
                     inject_kwargs={'bfgs_threshold': BFGSTHRESH})
    def __init__(self, *args, **kwargs):
        self.args = args
        self.kwargs = kwargs
        super(SmartMin, self).__init__(*args, **kwargs)
        self.bfgs_threshold = kwargs.get('bfgs_threshold', BFGSTHRESH)

    def run(self):
        descent_kwargs = self.kwargs.copy()
        descent_kwargs['force_tolerance'] = self.bfgs_threshold
        descender = GradientDescent(*self.args, **descent_kwargs)
        descender.run()

        if descender.current_step >= self.nsteps:
            self.traj = descender.traj
            return

        bfgs_kwargs = self.kwargs.copy()
        bfgs_kwargs['_restart_from'] = descender.current_step
        bfgs_kwargs['_restart_energy'] = descender._initial_energy
        bfgs_kwargs.setdefault('frame_interval', descender.frame_interval)
        bfgsmin = BFGS(*self.args, **bfgs_kwargs)
        bfgsmin.current_step = descender.current_step
        bfgsmin.run()

        self.traj = descender.traj + bfgsmin.traj
        self.current_step = bfgsmin.current_step


minimize = SmartMin._as_function('minimize')
exports(minimize)
toplevel(minimize)

