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
# TODO: look into http://molmod.github.io/

"""
"Generic" energy models - models that can be specified directly by name, without worrying about
which specific implementation is used.

Currently, everything here is an alias. However, more complicated logic (including runtime
dispatch) may be used to determine the best implementation in a given situation
"""
from moldesign import utils

from . import PySCFPotential
from . import OpenMMPotential


##################
# ForceField
@utils.exports
class ForceField(OpenMMPotential):
    pass  # currently an alias


##################
# QM generics
@utils.exports
class RHF(PySCFPotential):
    @utils.doc_inherit
    def __init__(self, *args, **kwargs):
        kwargs['theory'] = 'rhf'
        super(RHF, self).__init__(*args, **kwargs)


@utils.exports
class DFT(PySCFPotential):
    @utils.doc_inherit
    def __init__(self, *args, **kwargs):
        kwargs['theory'] = 'rks'
        super(DFT, self).__init__(*args, **kwargs)


@utils.exports
class B3LYP(PySCFPotential):
    @utils.doc_inherit
    def __init__(self, *args, **kwargs):
        kwargs['theory'] = 'rks'
        kwargs['functional'] = 'b3lyp'
        super(B3LYP, self).__init__(*args, **kwargs)


@utils.exports
class MP2(PySCFPotential):
    @utils.doc_inherit
    def __init__(self, *args, **kwargs):
        kwargs['theory'] = 'mp2'
        super(MP2, self).__init__(*args, **kwargs)

@utils.exports
class CASSCF(PySCFPotential):
    @utils.doc_inherit
    def __init__(self, *args, **kwargs):
        kwargs['theory'] = 'casscf'
        super(CASSCF, self).__init__(*args, **kwargs)