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
from __future__ import print_function, absolute_import, division

import os as _os

_building_docs = bool(_os.environ.get('SPHINX_IS_BUILDING_DOCS', ""))

from . import data
PACKAGEPATH = data.PACKAGEPATH

# Base subpackages - import these first
from . import utils
from . import units

# Functional subpackages
from . import compute
from . import fileio
from . import exceptions
from . import external
from . import forcefields
from . import geom
from . import helpers
from . import integrators
from . import interfaces
from . import parameters
from . import mathutils
from . import min
from . import models
from . import method
from . import orbitals
from . import molecules
from . import tools
from . import widgets

# Populate the top-level namespace (imports everything from each <submodule>.__all__ variable)
from .exceptions import *
from .fileio import *
from .forcefields import *
from .geom import *
from .min import *
from .orbitals import *
from .molecules import *
from .tools import *

# Set up cloud computing
compute.init_config()

# package metadata
from . import _version
__version__ = _version.get_versions()['version']
__copyright__ = "Copyright 2017 Autodesk Inc."
__license__ = "Apache 2.0"

# Set warnings appropriately
# TODO: don't clobber user's settings!!!
import numpy as _np
import warnings as _warnings
_np.seterr(all='raise')
_warnings.simplefilter('error', _np.ComplexWarning)

# For documentation purposes only - make sphinx document the toplevel namespace
if _building_docs:
    __all__ = fileio.__all__ + \
              geom.__all__ + \
              min.__all__ + \
              orbitals.__all__ + \
              molecules.__all__ + \
              tools.__all__

