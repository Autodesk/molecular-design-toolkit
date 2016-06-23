__all__ = []

def toplevel(o):
    __all__.append(o.__name__)
    return o

from .atoms import *
from .biounits import *
from .molecule import *
from .trajectory import *
from .converters import *

from moldesign import _building_docs
if not _building_docs:
    __all__ = (
               'read from_smiles from_pdb from_name build_bdna write read_amber '  # converters.py
               'build_assembly mol_to_openmm_sim mol_to_pybel '  # converters.py
               'Molecule NotCalculatedError NoConvergence MolecularProperties '  # molecule.py
               'Trajectory TrajectoryViewer '  # trajectory.py
               ).split()
