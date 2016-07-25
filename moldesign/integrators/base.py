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
from moldesign import parameters
from moldesign.method import Method


class IntegratorBase(Method):
    """Base class for all integrators"""

    PARAMETERS = parameters.integrator_parameters.values()

    def run(self, run_for):
        """
        To be called by parent molecule
        :param run_for: number of steps (if integer), or amount of time (if has units of time)
        :return: trajectory
        """
        raise NotImplementedError('This is an abstract base class!')

    def prep(self):
        """
        Prepare to run. Possibly do a test run.
        This might need to call the mol.model.build. Make sure you don't have a
        circular call here
        :return:
        """
        raise NotImplementedError()

    @staticmethod
    def time_to_steps(time, timestep):
        try:
            dims = time.dimensionality
            assert len(dims) == 1 and dims['[time]'] == 1.0
        except (AttributeError, AssertionError):
            assert type(time) == int, "argument to integrator.run must have units of time or be an int"
            return time
        else:
            return int(round(time / timestep))


class MDBase(IntegratorBase):
    PARAMETERS = IntegratorBase.PARAMETERS + parameters.md_parameters.values()


class ConstantTemperatureBase(MDBase):
    PARAMETERS = MDBase.PARAMETERS + parameters.constant_temp_parameters.values()


class LangevinBase(ConstantTemperatureBase):
    PARAMETERS = ConstantTemperatureBase.PARAMETERS + parameters.langevin_parameters.values()