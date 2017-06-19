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
from past.builtins import basestring
import future.utils

import os
import sys
import yaml
import warnings

from pyccc import engines

from .. import utils
from . import compute
from .. import _version

default_engine = None

RUNNING_ON_WORKER = (os.environ.get('IS_PYCCC_JOB', '0') == '1')

COMPUTE_CONNECTION_WARNING = """
WARNING: Failed to connect to a computational engine - MDT won't be able to run
anything outside of Python.

You can fix this either:
  1) interactively, in a notebook, by running `moldesign.configure()`, or,
  2) by modifying the configuration dictionary `moldesign.compute.config`, then
     running `moldesign.compute.reset_compute_engine()` to try again."""


config = utils.DotDict()
""" dict: dictionary of parameters (read from user's moldesign.yml at startup)

This dictionary configures how MDT interacts with "computational engines" - the systems that run
jobs, usually using docker images.


Notes:
    Values in this dictionary can be configured at runtime; however, after changing them,
    you should update the compute configuration by running
    ``moldesign.compute.reset_compute_engine()``

Configuration is specified using the following keys:

Args:
    engine_type (str): The computational job engine to use. Currently supported:
        'docker' (default) or 'subprocess'
    default_repository (str): Repository to pull MDT's standard images from.
        default: 'docker.io/autodesk/moldesign:'
    default_python_image (str): Image to run python commands in
        (default: ``docker.io/autodesk/moldesign:moldesign_complete-[VERSION]``, where
        [VERSION] is the version of MDT)
    default_version_tag (str): Default version tag for docker images
         (default: ``moldesign.__version__``)
    default_docker_url (str): URL for communicating with docker if ``engine_type=='docker'``.
        (default: Determined from $DOCKER_HOST, usually this will be the client you run on
        the command line)

MDT uses a non-standard docker tagging system to store its docker images. Generally,
a given image is pulled from a URL of the form:
   ``http://[some_docker_registry/orgname/reponame]:[imagename]-[versiontag]``

For instance, the ambertools 0.6 image is stored in the default repository at:
   ``http://docker.io/autodesk/moldesign:ambertools-0.6``
"""

ENVVAR = "MOLDESIGN_CONFIG"
""" str: name of environmental variable that stores path to moldesign.yml file.

If this variable is not set, ``$HOME/.moldesign/moldesign.yml`` will be used by default."""

DEFAULT_CONFIG_PATH = os.path.join(os.environ['HOME'], '.moldesign/moldesign.yml')
""" str: default search path for moldesign.yml."""

# TODO: we're currently hardcoding this at release - there's got to be a better way
DEFAULT_VERSION_TAG = _version.get_versions()['version']

CONFIG_DEFAULTS = utils.DotDict(engine_type='docker',
                                default_repository='docker.io/autodesk/moldesign:',
                                default_python_image=None,
                                default_docker_host='',
                                default_version_tag=DEFAULT_VERSION_TAG,
                                devmode=False)

DEF_CONFIG = CONFIG_DEFAULTS.copy()
""" dict: default configuration to be written to moldesign.yml if it doesn't exist
"""


def registry_login(client, login):
    print('Logging into docker registry @ %s ...' % login['registry'], end=' ')
    sys.stdout.flush()
    try:
        client.login(**login)
    except Exception as exc:
        warnings.warn('Failed to connect to %s. ' % login['registry'] +
                      'Some functions may time out.')
    else:
        print('done')


def write_config(path=None):
    if path is None:
        path = _get_config_path()

    if not os.path.exists(path) and path == DEFAULT_CONFIG_PATH:
        confdir = os.path.dirname(path)
        if not os.path.exists(confdir):
            os.mkdir(confdir)
            print('Created moldesign configuration directory %s' % confdir)

    configdump = {}
    for key in CONFIG_DEFAULTS:
        if key in config and isinstance(config[key], (basestring, int, float)):
            configdump[key] = config[key]
    with open(path, 'w') as f:
        yaml.safe_dump(configdump, f, default_flow_style=False)

    print('Wrote moldesign configuration to %s' % path)


def get_engine():
    from moldesign import compute
    if compute.default_engine is None:
        try:
            reset_compute_engine()
        except:
            print(COMPUTE_CONNECTION_WARNING)
            raise
    return compute.default_engine


def init_config():
    """Called at the end of package import to read initial configuration and setup cloud computing.
    """
    config.update(CONFIG_DEFAULTS)

    path = _get_config_path()
    if os.path.exists(path):
        try:
            with open(path, 'r') as infile:
                print('Reading configuration from %s' % path)
                newconf = yaml.load(infile)
                if not isinstance(newconf, dict):
                    raise TypeError('Cannot read configuration "%s" from %s.' % (newconf, path))
        except (IOError, KeyError, TypeError) as e:
            print(('WARNING: exception while reading configuration: %s. '
                   'using built-in default configuration') % e)
        else:
            config.update(newconf)

        _check_override('default_version_tag', DEFAULT_VERSION_TAG, path)

    if 'default_version_tag' not in config:
        config.default_version_tag = DEFAULT_VERSION_TAG

    if future.utils.PY2:
        expcted_docker_python_image = compute.get_image_path('moldesign_complete_py2')
    else:
        expcted_docker_python_image = compute.get_image_path('moldesign_complete')

    if config.get('default_python_image', None) is None:
        config.default_python_image = expcted_docker_python_image


def _check_override(tagname, expected, path):
    if tagname in config and config.default_version_tag != expected:
        print ('WARNING: Configuration file specifies a different value for %s! '
               "Remove the `%s` field from %s unless you know what you're doing"
               % (tagname, tagname, path))


def _get_config_path():
    if ENVVAR in os.environ:
        path = os.environ[ENVVAR]
    else:
        path = DEFAULT_CONFIG_PATH
    return path


def reset_compute_engine():
    """Read the configuration dict at ``moldesign.compute.config`` and set up compute engine

    Sets the module-level variable ``default_engine``.

    Returns:
        dict: copy of the config dictionary used to set the engine
    """
    from moldesign import compute

    compute.default_engine = None

    if config.engine_type == 'docker':
        if config.default_docker_host:
            notice = 'Connecting to docker host at %s' % config.default_docker_host
            hosturl = config.default_docker_host
        else:
            notice = "Connecting to your docker engine"
            hosturl = None
        with utils.textnotify(notice):
            compute.default_engine = engines.Docker(hosturl)
        _connect_docker_registry()

    elif config.engine_type == 'subprocess':
        compute.default_engine = engines.Subprocess()
        print("""WARNING: running all computational jobs as subprocesses on this machine.
This requires that you have all necessary software installed locally.
To change the engine, call moldesign.configure() or modify moldesign.compute.config .""")

    elif config.engine_type in ('ccc', 'docker-machine'):
        raise ValueError('Computational engine type "%s" is no longer supported by MDT. '
                         'Please install docker (https://docker.io) and set: \n'
                         '   engine_type: docker'
                         'in ~/.moldesign/moldesign.yml')

    else:
        raise ValueError('Unrecognized engine %s' % config.engine_type)


def _connect_docker_registry():
    from moldesign import compute

    if config.get('docker_registry_login', None):
        registry_login(compute.default_engine.client, config.docker_registry_login)
