from __future__ import absolute_import

import moldesign.__main__

from moldesign import data
PACKAGEPATH = data.PACKAGEPATH

from moldesign import basemethods, compute, data, helpers, integrators, interfaces, logs
from moldesign import minimizers, models, orbitals, symmetry, ui, utils

# Make objects from these modules accessible at the top level
# (note that they all use __all__ to specify what gets exported)
from moldesign import atoms, biounits, converters, forcefield, geometry
from moldesign import molecule, tools, trajectory

from moldesign.atoms import *
from moldesign.biounits import *
from moldesign.converters import *
from moldesign.forcefield import *
from moldesign.geometry import *
from moldesign.molecule import *
from moldesign.tools import *
from moldesign.trajectory import *

# This is here primarily for sphinx's benefit
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
