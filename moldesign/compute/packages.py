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
import subprocess
import importlib
import pkg_resources
from future.utils import PY2

from .. import utils


class InterfacedPackage(object):
    def __init__(self, packagename, expectedversion, docker_name='moldesign_complete',
                 importname=None, required=False):
        self.name = packagename  # for named_dict
        self.packagename = packagename
        if importname is None:
            self.importname = packagename
        else:
            self.importname = importname
        self.expectedversion = expectedversion
        self.required = required
        self.docker_name = docker_name

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

    @property
    def force_remote(self):
        from .configuration import config
        return config['run_remote'].get(self.name, False)

    @force_remote.setter
    def force_remote(self, value):
        from .configuration import config
        config['run_remote'][self.name] = value

    def runsremotely(self, _f=None, **kwargs):
        """ Wraps functions that can be run in a remote docker container.

        The function will run remotely whenever ``self.force_remote=True``
        """
        from .runsremotely import runsremotely
        from . import compute

        def should_run_remote():
            return self.force_remote

        def wrapper(f):
            wrapped = runsremotely(image=compute.get_image_path(self.docker_name),
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


class InterfacedExecutable(object):
    def __init__(self, exename, expectedversion, docker_image, path=None, version_flag='--version'):
        self.name = exename  # for named_dict
        self.exename = exename
        self.expectedversion = expectedversion
        self.image_name = docker_image
        self._path = path
        self.version_flag = version_flag

    @property
    def run_local(self):
        from .configuration import config
        return config['run_local'].get(self.name, False)

    @run_local.setter
    def run_local(self, value):
        from .configuration import config
        config['run_local'][self.name] = value

    @property
    def path(self):
        if self._path is not None:
            return self._path

        p = utils.which(self.exename)
        if p is not None:
            return p

        # check local directory too
        elif os.path.isfile(self.exename) and os.access(self.exename, os.X_OK):
            return os.path.abspath(self.exename)

        else:
            return None

    @path.setter
    def path(self, value):
        self._path = value

    def installed_version(self):
        if not self.is_installed():
            return None

        try:
            return subprocess.check_output([self.path, self.version_flag]).strip()
        except subprocess.CalledProcessError:
            return None

    def is_installed(self):
        if self.path is None:
            return False

        return os.path.exists(self.path)

    def make_job(self, **kwargs):
        import pyccc
        from . import compute

        kwargs['submit'] = False
        if self.run_local:
            kwargs['engine'] = pyccc.Subprocess()
        else:
            kwargs['image'] = compute.get_image_path(self.image_name)
        job = pyccc.Job(**kwargs)
        return job


nwchem = InterfacedExecutable('nwchem.exe', None, 'nwchem', version_flag=None)
opsin = InterfacedExecutable('opsin', None, 'opsin', version_flag=None)
nab = InterfacedExecutable('nab', '16', 'nucleic_acid_builder', version_flag=None)
symmol = InterfacedExecutable('symmol', None, 'symmol', version_flag=None)
tleap = InterfacedExecutable('tleap', '16', 'ambertools', version_flag=None)
antechamber = InterfacedExecutable('antechamber', '16', 'ambertools', version_flag=None)
nbo = InterfacedExecutable('nbo', '6.1', 'nbo', version_flag=None)

executables = [nwchem, opsin, nab, symmol, tleap, antechamber]


__all__ = [x.name for x in packages + executables]
