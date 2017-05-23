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

import moldesign as mdt

from .. import data
from ..utils import exports
from . import TLeapForcefield


@exports
class DefaultAmber(TLeapForcefield):
    def __init__(self):
        sourcelines = ['source %s' % data.AMBER_LEAPRC[ffname] for ffname in data.AMBER_DEFAULT]
        super(DefaultAmber, self).__init__(sourcelines)


@exports
class GAFF2ForceField(TLeapForcefield):
    def __init__(self):
        super(DefaultAmber, self).__init__(['source leaprc.gaff2'])
