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
# TODO: look into http://molmod.github.io/

"""
"Generic" energy models - models that can be specified directly by name, without worrying about
which specific implementation is used.

Currently, everything here is an alias. However, more complicated logic (including runtime
dispatch) may be used to determine the best implementation in a given situation
"""

from . import PySCFPotential
from . import OpenMMPotential


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


##################
# ForceField
@exports
def ForceField(**kwargs):
    return OpenMMPotential(**kwargs)


##################
# QM generics
@exports
def RHF(**kwargs):
    return PySCFPotential(theory='rhf', **kwargs)


@exports
def DFT(**kwargs):
    return PySCFPotential(theory='dft', **kwargs)


@exports
def B3LYP(**kwargs):
    return PySCFPotential(theory='dft', funtional='b3lyp', **kwargs)


@exports
def MP2(**kwargs):
    return PySCFPotential(theory='mp2', **kwargs)
