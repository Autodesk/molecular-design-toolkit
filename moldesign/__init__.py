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
import os as _os

_building_docs = bool(_os.environ.get('SPHINX_IS_BUILDING_DOCS', ""))

import moldesign.__main__

from . import data
PACKAGEPATH = data.PACKAGEPATH

# Import all subpackages / submodules
from . import compute
from . import converters
from . import exceptions
from . import external
from . import forcefields
from . import geom
from . import helpers
from . import integrators
from . import interfaces
from . import keywords
from . import mathutils
from . import min
from . import models
from . import method
from . import orbitals
from . import structure
from . import tools
from . import uibase
from . import units
from . import utils
from . import viewer
from . import widgets

# Populate the top-level namespace
from .converters import *
from .geom import *
from .orbitals import *
from .structure import *
from .viewer import *

# package metadata
from moldesign import _version
__version__ = _version.get_versions()['version']
__copyright__ = "Copyright 2016 Autodesk Inc."
__license__ = "Apache 2.0"

# This forces sphinx to document the top-level namespace
if _building_docs:
    __all__ = (atoms.__all__+
               biounits.__all__+
               converters.__all__+
               forcefield.__all__+
               coords.__all__+
               molecule.__all__+
               tools.__all__+
               trajectory.__all__)


# Set warnings appropriately
# TODO: don't clobber user's or other package's settings!!!
import numpy as _np
import warnings as _warnings
_np.seterr(all='raise')
_warnings.simplefilter('error', _np.ComplexWarning)


