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

import moldesign as mdt

from .. import data
from ..utils import exports
from . import TLeapForcefield


@exports
class TLeapLib(TLeapForcefield):
    """ A forcefield from TLeap's library
    """
    def __init__(self, *ffnames):
        self.names = ffnames
        sourcelines = ['source %s' % data.AMBER_LEAPRC[ffn] for ffn in ffnames]
        super().__init__(sourcelines)

    def __str__(self):
        return "TLeap data: %s" % ','.join(self.names)


@exports
class DefaultAmber(TLeapLib):
    def __init__(self):
        super().__init__(*data.AMBER_DEFAULT)


@exports
class GAFF2ForceField(TLeapLib):
    def __init__(self):
        super().__init__('gaff2')
