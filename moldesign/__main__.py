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
"""
This file collects the various command line tasks accessed via
``python -m moldesign [command]``

The functions here are intended help users set up their environment.

Note that MDT routines will NOT be importable from this module when it runs as a script -- you
won't be working with molecules or atoms in this module.
"""
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
EXAMPLE_DIR_SRC = unit_def_file = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                                               '_notebooks')


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

    elif args.command == 'config':
        print 'Reading config file from: %s' % CONFIG_PATH
        print '----------------------------'
        with open(CONFIG_PATH, 'r') as cfgfile:
            for key, value in yaml.load(cfgfile).iteritems():
                print '%s: %s' % (key, value)

    else:
        assert False

DOCKER_IMAGES = 'ambertools14 moldesign moldesign_notebook opsin symmol python_install'.split()
def pull():
    from moldesign import compute
    streams = []
    for img in DOCKER_IMAGES:
        imgurl = compute.get_image_path(img)
        print 'Pulling %s' % imgurl
        streams.append(compute.get_engine().client.pull(imgurl))
    for s in streams:
        for line in s.split('\n'):
            print line


def launch(cwd=None, path=''):
    server, portnum = launch_jupyter_server(cwd=cwd)
    wait_net_service('localhost', portnum)
    open_browser('http://localhost:%d/%s' % (portnum, path))
    try:
        server.wait()
    except KeyboardInterrupt:
        print 'Shutting down ...'
        server.terminate()
        server.wait()
        print 'Jupyter terminated.'


def copy_example_dir(use_existing=False):
    print 'Copying MDT examples to `%s` ...' % EXAMPLE_DIR_TARGET
    if os.path.exists(EXAMPLE_DIR_TARGET):
        if use_existing:
            return
        else:
            print '\n'.join(
                    ['FAILED: directory already exists. Please:'
                     ' 1) Rename or remove the existing directory at %s,' % EXAMPLE_DIR_TARGET,
                     ' 2) Run this command in a different location, or'
                     ' 3) Run `python -m moldesign intro` to launch the example gallery.'])
            sys.exit(1)
    else:
        shutil.copytree(EXAMPLE_DIR_SRC, EXAMPLE_DIR_TARGET)
        print 'Done.'


def launch_jupyter_server(cwd=None):
    for i in xrange(8888, 9999):
        if localhost_port_available(i):
            portnum = i
            break
    else:
        print 'WARNING: no available port found between 8888 and 9999. Will try a random port ... '
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
    print 'Point your browser to %s to get started.' % url  # fallback


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
            print 'x'
            # this exception occurs only if timeout is set
            if timeout: return False

        except socket.error as err:
            if type(err.args) != tuple or err[0] not in (errno.ETIMEDOUT, errno.ECONNREFUSED):
                raise
        else:
            s.close()
            return True

        print '.',
        sys.stdout.flush()
        time.sleep(0.1)


def check_path(exes):
    return {c: distutils.spawn.find_executable(c) for c in exes}

is_mac = (check_path(['osascript'])['osascript'] is not None)


if __name__ == '__main__':
    main()
