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
from moldesign.utils import exports

from .base import MMBase
from .jsonmodel import JsonModelBase

IMAGE = 'nwchem'


@exports
class SanderMM(JsonModelBase, MMBase):
    """ Interface with NWChem package (QM only)

    Note:
        This is the first interface based on our new wrapping strategy. This is slightly hacked,
        but has the potential to become very general; very few things here are NWChem-specific
    """
    IMAGE = 'ambertools'
    MODELNAME = 'sander'
    DEFAULT_PROPERTIES = ['potential_energy']
    ALL_PROPERTIES = DEFAULT_PROPERTIES + ['forces']

    RUNNER = 'runsander.py'
    PARSER = 'parsesander.py'
