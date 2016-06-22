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
import os
import sys

import yaml
import warnings
from pyccc import engines, job
import pyccc.python as pycccpy

from moldesign.utils.classes import DotDict

cfg = launchcmd = launchpy = config_yaml = _docker_client = None
ENVVAR = "MOLDESIGN_CONFIG"
DEFAULT_CONFIG_PATH = os.path.join(os.environ['HOME'], '.moldesign/moldesign.yml')
DEF_CONFIG = DotDict(default_engine='docker',
                     default_repository='docker-hub.autodesk.com/virshua/moldesign:',
                     default_image='moldesign_complete',
                     default_docker_url='unix://var/run/docker.sock',
                     default_docker_machine=None,
                     default_num_threads=4,
                     version_tag='latest',
                     show_splash=True,
                     _splash_done=False)


def registry_login(client, login):
    print 'Logging into docker registry @ %s ...' % login['registry'],
    sys.stdout.flush()
    try:
        client.login(**login)
    except Exception as exc:
        warnings.warn('Failed to connect to %s. ' % login['registry'] +
                      'Some functions may time out.')
    else:
        print 'done'


def print_configuration():
    print 'Computing platform:', cfg.default_engine
    print 'Platform URL:', cfg.get('default_docker_url', 'none')
    print 'Image location:', cfg.default_repository
    print 'Image tag:', cfg.version_tag
    print 'Image for python commands:', cfg.default_image


def _get_config():
    # TODO: refactor, allow runtime_ms changes, move host setup to PyCCC. see BBLL-64
    # TODO: how to configure workers? (especially ones that don't need to launch any jobs themselves)
    """Called at import to configure this module"""
    global cfg, launchcmd, launchpy, config_yaml, _docker_client
    config_yaml = DotDict(DEF_CONFIG)
    cfg = DotDict(config_yaml)

    if ENVVAR in os.environ:
        path = os.environ[ENVVAR]
    else:
        path = DEFAULT_CONFIG_PATH

    if os.path.exists(path):
        with open(path, 'r') as infile:
            print 'Reading configuration from %s' % path
            config_yaml.update(yaml.load(infile))
    else:
        print 'No config file found at %s - using defaults' % path

    if config_yaml.default_engine == 'docker':
        import docker
        from pyccc import docker_utils
        if config_yaml.get('default_docker_machine', None) is not None:
            client = docker_utils.docker_machine_client(config_yaml.default_docker_machine)
            cfg.default_docker_url = client.base_url
        else:
            client = docker.Client(base_url=config_yaml.default_docker_url)
        cfg.default_engine = engines.Docker(client=client)
        _docker_client = client

        if 'docker_registry_login' in config_yaml:
            registry_login(client, config_yaml.docker_registry_login)
    elif config_yaml.default_engine == 'subprocess':
        cfg.default_engine = engines.Subprocess
    else:
        raise ValueError('Unrecognized engine %s' % config_yaml.default_engine)

    cfg.default_image = config_yaml.default_image
    cfg.default_repository = config_yaml.default_repository
    cfg.version_tag = config_yaml.version_tag
    launchcmd = job.Launcher(cfg.default_engine,
                             cfg.default_image)
    launchpy = pycccpy.PythonLauncher(cfg.default_engine,
                                      cfg.default_image)


_get_config()
