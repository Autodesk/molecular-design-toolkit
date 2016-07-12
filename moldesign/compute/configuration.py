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
import subprocess
import yaml
import warnings

from pyccc import engines

from moldesign import utils

default_engine = None

FREE_COMPUTE_CANNON = 'cloudcomputecannon.bionano.autodesk.com:9000'

RUNNING_ON_WORKER = (os.environ.get('IS_PYCCC_JOB', '0') == '1')


# TODO: *write* configuration plus initial install default; sensible defaults

config = utils.DotDict()
""" dict: dictionary of parameters (read from user's moldesign.yml at startup)

This dictionary configures how MDT interacts with "computational engines" - the systems that run
jobs, usually using docker images.


Notes:
    Values in this dictionary can be configured at runtime; however, after changing them,
    you should update the compute configuration by running ``moldesign.compute.reset_from_config()``

Configuration is specified using the following keys:

Args:
    engine_type (str): The computational job engine to use. Default: 'ccc'. Currently supported:
        'docker', 'ccc', or 'subprocess'
    default_repository (str): Repository to pull MDT's standard images from.
        default: 'docker.io/autodesk/moldesign:'
    default_python_image (str): Image to run python commands in
        (default: ``docker.io/autodesk/moldesign:moldesign_complete-[VERSION]``, where
        [VERSION] is the version of MDT)
    default_version_tag (str): Default version tag for docker images (default: 'latest')
    default_docker_machine (str): Name of the docker machine to connect to; if
        ``engine_type=='docker'``, EITHER this OR ``default_docker_url`` (but not both) must
        be set. (default: 'default')
    default_docker_url (str): URL for the docker daemon to run; if ``engine_type=='docker'``,
        EITHER this OR ``default_docker_machine`` (but not both) must be set.
        (default: ``unix://var/run/docker.sock``, the URL for a local docker instance)

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


CONFIG_DEFAULTS = utils.DotDict(engine_type='docker-machine',
                                default_repository='docker-hub.autodesk.com/virshua/moldesign:',
                                default_ccc_host=FREE_COMPUTE_CANNON,
                                default_python_image='moldesign_complete',
                                default_docker_host='unix://var/run/docker.sock',
                                default_docker_machine='default',
                                default_version_tag='latest')

DEF_CONFIG = CONFIG_DEFAULTS.copy()
""" dict: default configuration to be written to moldesign.yml if it doesn't exist
"""
for x in 'default_docker_machine'.split(): DEF_CONFIG.pop(x)


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


def init_config():
    """Called at the end of package import to read initial configuration and setup cloud computing.

    At runtime, call :meth:`reset_from_config` to change the configuration
    """
    config.update(CONFIG_DEFAULTS)

    if ENVVAR in os.environ:
        path = os.environ[ENVVAR]
    else:
        path = DEFAULT_CONFIG_PATH

    if os.path.exists(path):
        with open(path, 'r') as infile:
            print 'Reading configuration from %s' % path
            config.update(yaml.load(infile))
    else:
        print 'No config file found at %s - using defaults' % path

    reset_from_config()


def reset_from_config():
    """Read the configuration dict at ``moldesign.compute.config`` and set up compute engine

    Returns:
        dict: copy of the config dictionary (for posterity)
    """
    from moldesign import compute

    if config.engine_type in ('docker', 'docker-machine'):
        compute.default_engine = _setup_docker()

    elif config.engine_type == 'subprocess':
        compute.default_engine = engines.Subprocess

    else:
        raise ValueError('Unrecognized engine %s' % config.default_engine)


def _setup_docker():
    """ Connect to a docker server, either directly or via docker-machine
    """
    import docker
    from pyccc import docker_utils
    if config.engine_type == 'docker-machine':
        dm = config.default_docker_machine
        success = False
        try:
            status = subprocess.check_output(['docker-machine', 'status', dm]).strip()
        except (subprocess.CalledProcessError, OSError):
            print 'WARNING: could not find docker-machine named "%s"' % dm
        else:
            if status != 'Running':
                print 'WARNING: docker-machine %s returned status: "%s"' % (dm, status)
            else:
                try:
                    with utils.textnotify('Connecting to docker-machine "%s"' % dm):
                        dockerclient = docker_utils.docker_machine_client(dm)
                except (subprocess.CalledProcessError, OSError) as e:
                    print 'WARNING: error connecting to docker-machine "%s": %s' % (dm, e)
                else:
                    success = True

        if not success:
            print 'WARNING: failed to connect to docker machine - cannot run jobs.'
            print 'After changing the configuration (moldesign.config.default_docker_machine) ' \
                  'and/or starting the docker-machine, run ' \
                  'moldesign.compute.reset_from_config() to try again.'
            return None

    elif config.engine_type == 'docker':
        with utils.textnotify('Connecting to docker at ' % config.default_docker_url):
            dockerclient = docker.Client(base_url=config.default_docker_url)

    if 'docker_registry_login' in config:
        registry_login(dockerclient, config.docker_registry_login)

    return engines.Docker(client=dockerclient)
