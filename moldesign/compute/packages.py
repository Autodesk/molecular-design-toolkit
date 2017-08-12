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
from future.utils import PY2
import importlib
import pkg_resources


class InterfacedPackage(object):
    def __init__(self, packagename, expectedversion, docker_image='moldesign_complete',
                 importname=None, required=False):
        self.name = packagename  # for named_dict
        self.packagename = packagename
        if importname is None:
            self.importname = packagename
        else:
            self.importname = importname
        self.expectedversion = expectedversion
        self.required = required
        self.force_remote = not self.is_installed()  # user can change this manually later
        self.docker_image = docker_image

    if PY2:
        def is_installed(self):
            import imp
            try:
                imp.find_module(self.importname)
            except (ImportError, OSError) as exc:
                return False
            else:
                return True

    else:
        def is_installed(self):
            spec = importlib.util.find_spec(self.importname)
            return spec is not None

    def installed_version(self):
        try:
            return pkg_resources.get_distribution(self.packagename).version
        except pkg_resources.DistributionNotFound:
            return None

    def runsremotely(self, _f=None, **kwargs):
        """ Wraps functions that can be run in a remote docker container.

        The function will run remotely whenever ``self.force_remote=True``
        """
        from .runsremotely import runsremotely

        def should_run_remote():
            return self.force_remote

        def wrapper(f):
            wrapped = runsremotely(image=self.docker_image,
                                   should_run_remote=should_run_remote,
                                   **kwargs)(f)
            return wrapped

        if _f is not None:
            assert not kwargs
            return wrapper(_f)
        else:
            return wrapper


biopython = InterfacedPackage('biopython', '1.68', importname='Bio', required=True)
parmed = InterfacedPackage('parmed', '2.7.3', required=True)

# can't find any run-time mechanism to get the version for openbabel ...
openbabel = InterfacedPackage('openbabel', '2.4')
pdbfixer = InterfacedPackage('pdbfixer', '1.4')
pyscf = InterfacedPackage('pyscf', '1.1')
openmm = InterfacedPackage('OpenMM', '7.1.1', importname='simtk')

packages = [biopython, parmed, openbabel, pdbfixer, pyscf, openmm]

__all__ = [x.name for x in packages]
