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
import os
import json

import moldesign as mdt
from .. import units as u

PACKAGEPATH = os.path.abspath(os.path.dirname(mdt.__file__))


class CyclicList(list):
    def __getitem__(self, item):
        return super().__getitem__(item % len(self))


COLOR_LIST = ['lightgreen', 'lightblue', 'lightgrey',
              'yellow', 'orange', 'purple', 'IndianRed',
              'PaleTurquoise', 'OldLace', 'Thistle', 'pink']


DEFAULT_FORCE_TOLERANCE = (0.0001 * u.hartree / u.bohr).defunits()  # taken from GAMESS OPTTOL keyword


def print_environment():
    """For reporting bugs - spits out the user's environment"""
    import sys
    version = {}
    for pkg in 'moldesign IPython ipywidgets jupyter matplotlib numpy docker pyccc distutils' \
               'nbmolviz jupyter_client jupyter_core pint Bio openbabel simtk pyscf pip setuptools'\
            .split():
        try:
            module = __import__(pkg)
        except ImportError as e:
            version[pkg] = str(e)
        else:
            try:
                version[pkg] = module.__version__
            except AttributeError as e:
                version[pkg] = str(e)
    env = {'platform': sys.platform,
           'version': sys.version,
           'prefix': sys.prefix}

    try:
        import platform
        env['machine'] = platform.machine()
        env['linux'] = platform.linux_distribution()
        env['mac'] = platform.mac_ver()
        env['windows'] = platform.win32_ver()
        env['impl'] = platform.python_implementation()
        env['arch'] = platform.architecture()
        env['system'] = platform.system()
        env['python_build'] = platform.python_build()
        env['platform_version'] = platform.version()

    except Exception as e:
        env['platform_exception'] = str(e)

    print(json.dumps({'env': env,
                      'versions': version}))
