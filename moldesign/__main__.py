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
"""
This file collects the various command line tasks accessed via
``python -m moldesign [command]``

The functions here are intended help users set up their environment.

Note that MDT routines will NOT be importable from this module when it runs as a script -- you
won't be working with molecules or atoms in this module.
"""
from __future__ import print_function
from builtins import range
import argparse
import distutils.spawn
import errno
import os
import random
import shutil
import socket
import subprocess
import sys
import time

import yaml

URL_OPENERS = ['Open', 'xdg-open', 'sensible-browser', 'gnome-open', 'x-www-browser']
JUPYTERPORT = 8888
DOCKER_TOOLS = 'docker docker-machine docker-compose'.split()
DOCKER_REPOSITORY = 'docker-hub.autodesk.com/virshua/moldesign:'
HOME = os.environ['HOME']
CONFIG_DIR = os.path.join(HOME, '.moldesign')
EXAMPLE_DIR_TARGET = os.path.join(os.path.curdir, 'moldesign-examples')
MOLDESIGN_SRC = os.path.abspath(os.path.dirname(__file__))
EXAMPLE_DIR_SRC = unit_def_file = os.path.join(MOLDESIGN_SRC, '_notebooks')
MDTVERSION = subprocess.check_output(['python', '-c',
                                      "import _version; print(_version.get_versions()['version'])"],
                                     cwd=MOLDESIGN_SRC).splitlines()[-1].decode('ascii')
VERFILEPATH = os.path.join(EXAMPLE_DIR_TARGET, '.mdtversion')


APPLESCRIPT_INSTALL_DOCKER = ('set response to (display dialog '
                              '"Molecular Design Toolkit needs the Docker Toolbox. Would you like '
                              'to download it now?" '
                              'buttons {"Quit", "Open download page"} '
                              'default button 2 with title "Buckyball Setup")\n'
                              'if button returned of response = "Open Download Page" then '
                              'open location "https://www.docker.com/docker-toolbox"'
                              'end if')
CONFIG_PATH = os.path.join(CONFIG_DIR, 'moldesign.yml')


def main():
    print('Molecular Design Toolkit v%s Launcher' % MDTVERSION)

    global CONFIG_PATH
    parser = argparse.ArgumentParser('python -m moldesign')

    subparsers = parser.add_subparsers(title='command', dest='command')
    subparsers.add_parser('intro', help='copy examples into current directory and launch a '
                                        'notebook')
    subparsers.add_parser('launch', help='launch a notebook and open it in a browser '
                                         '(equivalent to running "jupyter notebook")')
    subparsers.add_parser('pull', help='download docker containers that MDT requires ('
                                       'only when a docker client is configured)')
    subparsers.add_parser('config', help='print configuration and exit')
    subparsers.add_parser('copyexamples', help='Copy example notebooks')
    subparsers.add_parser('devbuild', help='rebuild required docker containers locally')
    subparsers.add_parser('devpull', help='Pull development images for latest release')
    subparsers.add_parser('version', help='Write version string and exit')

    parser.add_argument('-f', '--config-file', type=str,
                        help='Path to config file')

    args = parser.parse_args()

    if args.config_file:
        CONFIG_PATH = args.config_file

    if args.command == 'intro':
        copy_example_dir(use_existing=True)
        launch(cwd=EXAMPLE_DIR_TARGET,
               path='notebooks/Getting%20Started.ipynb')

    elif args.command == 'pull':
        pull()

    elif args.command == 'launch':
        launch()

    elif args.command == 'devbuild':
        devbuild()

    elif args.command == 'devpull':
        devpull()

    elif args.command == 'copyexamples':
        copy_example_dir(use_existing=False)

    elif args.command == 'version':
        print(MDTVERSION)

    elif args.command == 'config':
        print('Reading config file from: %s' % CONFIG_PATH)
        print('----------------------------')
        with open(CONFIG_PATH, 'r') as cfgfile:
            for key, value in yaml.load(cfgfile).items():
                print('%s: %s' % (key, value))

    else:
        raise ValueError("Unhandled CLI command '%s'" % args.command)

DOCKER_IMAGES = 'ambertools moldesign_complete opsin symmol nwchem'.split()

def pull():
    for img in DOCKER_IMAGES:
        _pull_img(img)


def _pull_img(img):
    from moldesign import compute
    imgurl = compute.get_image_path(img, _devmode=False)
    print('Pulling %s' % imgurl)
    subprocess.check_call(['docker', 'pull', imgurl])


BUILD_FILES = "nwchem_build pyscf_build".split()

def devbuild():
    print('-' * 80)
    print("Molecular design toolkit is downloading and building local, up-to-date copies of ")
    print("all docker containers it depends on. To use them, set:")
    print("   devmode: true")
    print("in ~/.moldesign/moldesign.yml.")
    print(('-' * 80) + '\n')

    devpull()

    subprocess.check_call('docker-make --all --tag dev'.split(),
                          cwd=os.path.join(MOLDESIGN_SRC, '..', 'DockerMakefiles'))


def devpull():
    for img in DOCKER_IMAGES+BUILD_FILES:
        try:
            _pull_img(img)
        except subprocess.CalledProcessError:
            print('"%s " not found. Will rebuild locally ...'%img)


def launch(cwd=None, path=''):
    server, portnum = launch_jupyter_server(cwd=cwd)
    wait_net_service('localhost', portnum)
    open_browser('http://localhost:%d/%s' % (portnum, path))
    try:
        server.wait()
    except KeyboardInterrupt:
        print('Shutting down ...')
        server.terminate()
        server.wait()
        print('Jupyter terminated.')


def copy_example_dir(use_existing=False):
    print('Copying MDT examples to `%s` ...' % EXAMPLE_DIR_TARGET)
    if os.path.exists(EXAMPLE_DIR_TARGET):
        check_existing_examples(use_existing)
    else:
        shutil.copytree(EXAMPLE_DIR_SRC, EXAMPLE_DIR_TARGET)
        with open(VERFILEPATH, 'w') as verfile:
            print(MDTVERSION, file=verfile)
        print('Done.')


def check_existing_examples(use_existing):
    if os.path.exists(VERFILEPATH):
        with open(VERFILEPATH, 'r') as vfile:
            version = vfile.read().strip()
    else:
        version = 'pre-0.7.4'

    if version != MDTVERSION:
        print ('WARNING - your example directory is out of date! It corresponds to MDT version '
               '%s, but you are using version %s'%(version, MDTVERSION))
        print ('If you want to update your examples, please rename or remove "%s"'
               % EXAMPLE_DIR_TARGET)

    if use_existing:
        return
    else:
        print('\n'.join(
                ['FAILED: directory already exists. Please:'
                 ' 1) Rename or remove the existing directory at %s,'%EXAMPLE_DIR_TARGET,
                 ' 2) Run this command in a different location, or'
                 ' 3) Run `python -m moldesign intro` to launch the example gallery.']))
        sys.exit(1)


def launch_jupyter_server(cwd=None):
    for i in range(8888, 9999):
        if localhost_port_available(i):
            portnum = i
            break
    else:
        print('WARNING: no available port found between 8888 and 9999. Will try a random port ... ')
        portnum = random.randint(10001, 65535)

    server = subprocess.Popen(('jupyter notebook --no-browser --port %d' % portnum).split(),
                              cwd=cwd)

    return server, portnum


def open_browser(url):
    for exe in URL_OPENERS:
        if distutils.spawn.find_executable(exe) is not None:
            try:
                subprocess.check_call([exe, url])
            except subprocess.CalledProcessError:
                continue
            else:
                return
    print('Point your browser to %s to get started.' % url)  # fallback


def localhost_port_available(portnum):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.settimeout(0.2)
    try:
        s.connect(("localhost", portnum))
    except socket.error as err:
        if err.errno == errno.ECONNREFUSED:
            return True
        else:
            raise
    else:
        return False


def yaml_dumper(*args):
    return yaml.dump(*args, default_flow_style=False)


def wait_net_service(server, port, timeout=None):
    """
    Wait for network service to appear
    FROM http://code.activestate.com/recipes/576655-wait-for-network-service-to-appear/
    @param timeout: in seconds, if None or 0 wait forever
    @return: True of False, if timeout is None may return only True or
             throw unhandled network exception
    """
    import socket
    import errno

    if timeout:
        from time import time as now
        # time module is needed to calc timeout shared between two exceptions
        end = now() + timeout

    while True:
        s = socket.socket()

        try:
            if timeout:
                next_timeout = end - now()
                if next_timeout < 0:
                    return False
                else:
                    s.settimeout(next_timeout)
            s.connect((server, port))

        except socket.timeout as err:
            print('x')
            # this exception occurs only if timeout is set
            if timeout: return False

        except socket.error as err:
            if type(err.args) != tuple or err[0] not in (errno.ETIMEDOUT, errno.ECONNREFUSED):
                raise
        else:
            s.close()
            return True

        print('.', end=' ')
        sys.stdout.flush()
        time.sleep(0.1)


def check_path(exes):
    return {c: distutils.spawn.find_executable(c) for c in exes}

is_mac = (check_path(['osascript'])['osascript'] is not None)


if __name__ == '__main__':
    main()
