def toplevel(o):
    __all__.append(o.__name__)
    return o
__all__ = []

from . import base
from .bfgs import *
from .descent import *
from .smart import *
