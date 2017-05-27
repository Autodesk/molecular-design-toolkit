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
import inspect

__all__ = 'exports exports_names'.split()


def exports(f):
    """ Add a function to its module's __all__ attribute
    """
    all_list = _get_module_all(1)
    all_list.append(f.__name__)
    return f


def exports_names(*names):
    """ Add names to this module's __all__ attribute
    """
    all_list = _get_module_all(1)
    all_list.extend(names)


def _get_module_all(depth):
    """
    Get a reference to the __all__ attribute of the module calling the function that
    calls this function :)
    FROM http://stackoverflow.com/q/6187355/1958900"""
    frm = inspect.stack()[depth+1]
    mod = inspect.getmodule(frm[0])
    if not hasattr(mod, '__all__'):
        mod.__all__ = []
    return mod.__all__

