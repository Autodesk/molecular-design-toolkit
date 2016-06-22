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
from __future__ import absolute_import

import os as _os

_building_docs = bool(_os.environ.get('SPHINX_IS_BUILDING_DOCS', ""))

import moldesign.__main__

from moldesign import data
PACKAGEPATH = data.PACKAGEPATH

from moldesign.widgets import logs
from moldesign.methods import models, forcefield, basemethods, integrators

# expose these modules at the top level -- TODO: don't?
from moldesign.structure import atoms, biounits, molecule, trajectory, converters
from moldesign import biounits, units, tools, data
from moldesign import converters
from moldesign.methods import forcefield
from moldesign.geom import geometry
from moldesign.structure import molecule
from moldesign import trajectory

# (note that they all use __all__ to specify what gets exported)
from moldesign.structure.atoms import *
from moldesign.structure.biounits import *
from moldesign.structure.converters import *
from moldesign.methods.forcefield import *
from moldesign.geom.geometry import *
from moldesign.structure.molecule import *
from moldesign.tools import *
from moldesign.structure.trajectory import *

# package metadata
from moldesign import _version
__version__ = _version.get_versions()['version']
__copyright__ = "Copyright 2016 Autodesk Inc."
__license__ = "Apache 2.0"

# This forces sphinx to document the top-level namespace
if _building_docs:
    __all__ = (atoms.__all__ +
               biounits.__all__ +
               converters.__all__ +
               forcefield.__all__ +
               geometry.__all__ +
               molecule.__all__ +
               tools.__all__ +
               trajectory.__all__)

# Set warnings appropriately
# TODO: don't clobber user's or other package's settings!!!
import numpy as _np
import warnings as _warnings
_np.seterr(all='raise')
_warnings.simplefilter('error', _np.ComplexWarning)


