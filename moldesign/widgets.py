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
import sys
import textwrap
import functools
from .utils import exports, exports_names

__all__ = 'BondSelector GeometryBuilder ResidueSelector Symmetrizer AtomSelector'.split()


try:
    import nbmolviz.widget_utils
except ImportError:
    nbmolviz_enabled = False
    nbmolviz_installed = False
else:
    nbmolviz_installed = True
    nbmolviz_enabled = nbmolviz.widget_utils.can_use_widgets()


def notebook_only_method(*args, **kwargs):
    assert not nbmolviz_enabled
    if nbmolviz_installed:
        raise ImportError("This function is only available in a Jupyter notebook!")
    else:
        raise ImportError(
                "The `nbmolviz` library must be installed to use this function!\n"
                "To use MDT's notebook visualization features, please install\n"
                "the nbmolviz library by running these commands at the terminal:\n"
                "  1) pip install nbmolviz\n"
                "  2) python -m nbmolviz activate\n"
                "Afterwards, restart and reload any running notebooks.")


if nbmolviz_enabled:
    # TODO: Make these into lazy imports
    from nbmolviz.widgets import (BondSelector, GeometryBuilder, ResidueSelector, Symmetrizer,
                                  AtomSelector)
    from nbmolviz.mdtconfig.compute import configure
    about = configure
    nbmolviz.widget_utils.print_extension_warnings(stream=sys.stderr)
else:
    BondSelector = GeometryBuilder = ResidueSelector = AtomSelector = Symmetrizer = configure \
        = about = notebook_only_method


def _get_nbmethod(name):
    # don't import nbmolviz methods until a method is actually called
    from nbmolviz import methods as nbmethods
    module = nbmethods
    for item in name.split('.'):
        module = getattr(module, item)
    return module


@exports
class WidgetMethod(object):
    def __init__(self, name):
        self.name = name
        self.method = None

    def __get__(self, instance, owner):
        if not nbmolviz_enabled:
            return notebook_only_method
        elif self.method is None:
            self.method = _get_nbmethod(self.name)
        return functools.partial(self.method, instance)

