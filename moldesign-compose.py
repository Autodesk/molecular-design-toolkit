#!/usr/bin/env python
"""
Standalone executable to deploy and run the latest version of moldesign.
"""
import atexit
import distutils.spawn
import argparse
import shutil
import os
import subprocess
import sys
import time

try:
    import yaml
    def yaml_dumper(*args):
        return yaml.dump(*args, default_flow_style=False)
except ImportError:
    import json
    yaml_dumper = json.dump


def main():
    if args.command == 'up':
        up()
    elif args.command in ['create', 'build', 'pull']:
        create()
    elif args.command in ['down', 'stop', 'kill']:
        down(args.command)
    elif args.command == 'rm':
        rm()
    elif args.command == 'status':
        print subprocess.check_output('docker-machine status %s' % args.docker_machine_name,
                                      shell=True).strip()
    else:
        raise ValueError('Unknown command %s' % args.command)    


def up():
    create()
    url = start_docker_machine(args.docker_machine_name)
    pull_image()
    process = compose_up(url, args.docker_machine_name)
    open_jupyter(url)
    try:
        process.wait()
    except KeyboardInterrupt:
        print 'Shutting down ...'


def create():
    check_docker_tools()
    make_config_dir(overwrite=args.delete_config)
    create_docker_machine(args.docker_machine_name)


def config():
    make_config_dir(overwrite=True)


def down(cmd):
    if machine_is_running():
        compose_down(cmd)
        subprocess.check_output('docker-machine stop %s' % args.docker_machine_name, shell=True)


def rm():
    down('kill')
    subprocess.check_call('docker-machine rm %s' % args.docker_machine_name, shell=True)
    

HOME = os.environ['HOME']
CONFIG_DIR = os.path.join(HOME, '.moldesign')
JUPYTERPORT = 8899
DOCKER_TOOLS = 'docker docker-machine docker-compose'.split()
DOCKER_COMPOSE_PATH = os.path.join(CONFIG_DIR, 'docker-compose.yml')
HOME = os.environ['HOME']


APPLESCRIPT_INSTALL_DOCKER = ('set response to (display dialog '
                              '"Buckyball needs the Docker Toolbox. Would you like to download it now?" '
                              'buttons {"Quit", "Open download page"} '
                              'default button 2 with title "Buckyball Setup")\n'
                              'if button returned of response = "Open Download Page" then '
                              'open location "https://www.docker.com/docker-toolbox"'
                              'end if')
CONFIG_PATH = os.path.join(CONFIG_DIR, 'moldesign.yml')


def machine_is_running():
    result = subprocess.check_output('docker-machine status %s' % args.docker_machine_name, shell=True)
    return result.strip() == 'Running'


def make_parser():
    p = argparse.ArgumentParser('Buckyball')
    p.add_argument('command', choices='up down rm config kill build pull restart status'.split())
    p.add_argument('--delete-config', action='store_true',
                        help='Write fresh configuration files to ~/.moldesign')
    p.add_argument('--version-tag',
                   help='docker registry version',
                   default='stable')
    p.add_argument('--repository',
                   help='Docker repository',
                   default='docker-hub.autodesk.com/virshua/moldesign:')
    
    vmargs = p.add_argument_group('VM options')
    vmargs.add_argument('-n', '--docker-machine-name', default='moldesign-app')
    vmargs.add_argument('-m', '--docker-machine-memory', default='8196',
                        help='Amount of RAM to allocate to newly created docker machine')
    vmargs.add_argument('-s', '--docker-machine-disk-size', default='50000',
                        help='Max VM disk size (MB)')
    vmargs.add_argument('-c', '--numcpus', default=-1,
                        help='Number of CPUs available to the VM (-1 implies all CPUs)')
    vmargs.add_argument('-d', '--storage-driver', default='overlay',
                        help='Docker storage backend driver (default: overlay)')
    vmargs.add_argument('-p', '--local-process-cap', default=8, type=int)
    return p


def check_path(exes):
    return {c: distutils.spawn.find_executable(c) for c in exes}


def image_name(base):
    if args.repository[-1] == ':':
        fmt = "{repo}{base}-{version}"
    elif args.repository[-1] == '/':
        fmt = "{repo}{base}:{version}"
    else:
        fmt = "{repo}/{base}:{version}"

    return fmt.format(repo=args.repository, version=args.version_tag, base=base)

def make_config_dir(overwrite=False):
    composefile = {
        "mdt": {"image": image_name('moldesign_notebook'),
               "ports":
                   ["8899:8888"],
               "environment":
                   {"BUCKYBALL_CONFIG": "%s/.moldesign/moldesign.yml" % HOME},
               "volumes": ["/var/run/docker.sock:/var/run/docker.sock",
                           "/Users:/Users",
                           "%s/mdtnotebooks:/notebooks/saved" % HOME]}}

    configfile = dict(default_provider='docker',
                      default_repository=args.repository,
                      default_docker_url='unix://var/run/docker.sock',
                      default_image=image_name("moldesign"),
                      version_tag=args.version_tag)

    if not os.path.exists(CONFIG_DIR):
        print 'Creating %s' % CONFIG_DIR
        os.makedirs(CONFIG_DIR)

    if overwrite or (not os.path.exists(DOCKER_COMPOSE_PATH)):
        print 'Creating default %s' % DOCKER_COMPOSE_PATH
        with open(DOCKER_COMPOSE_PATH, 'w') as outfile:
            yaml_dumper(composefile, outfile)

    if overwrite or (not os.path.exists(CONFIG_PATH)):
        print 'Creating default %s' % CONFIG_PATH
        with open(CONFIG_PATH, 'w') as outfile:
            yaml_dumper(configfile, outfile)


def check_docker_tools():
    exe_paths = check_path(DOCKER_TOOLS)
    missing_docker = [d for d in exe_paths if exe_paths[d] is None]
    if len(missing_docker) == 0:
        versions = {d: subprocess.check_output([d, '--version']).strip() for d in exe_paths}
        return versions

    # If here, prompt the user to install docker toolbox
    if is_mac:
        subprocess.check_call(['osascript', '-e', APPLESCRIPT_INSTALL_DOCKER])
    else:
        print 'Buckyball needs the Docker Toolbox. Please download and install it' \
              'from https://www.docker.com/docker-toolbox'
    sys.exit(1)


def start_docker_machine(name):
    status = subprocess.check_output(['docker-machine', 'status', name]).strip()
    if status != 'Running':
        print 'Starting docker machine:'
        cmd = 'docker-machine start {0}'.format(name)
        print cmd
        subprocess.check_call(cmd.split())
    else:
        print '%s is already running' % name

    set_machine_env(name)
    url = subprocess.check_output(['docker-machine', 'ip', name]).strip()
    return url


def set_machine_env(name):
    result = subprocess.check_output('docker-machine env %s' % name, shell=True)
    for line in result.split('\n'):
        fields = line.split()
        if len(fields) == 0: continue
        if fields[0] == 'export':
            key, val = fields[1].split('=')
            os.environ[key] = _dequote(val)  # strips the quotes from the string


def _dequote(s):
    if (s[0] == s[-1]) and s.startswith(("'", '"')):
        return s[1:-1]
    return s


def create_docker_machine(name):
    print '\n_____________________________________\nChecking for docker machine "%s"' % name
    try:
        subprocess.check_call(['docker-machine', 'status', name])
    except subprocess.CalledProcessError:
        pass
    else:
        print '%s already exists' % name
        return
    
    print '\n_____________________________________\nCreating docker machine: '
    cmd = ('docker-machine create -d virtualbox '
           '--engine-storage-driver {fs} '
           '--virtualbox-cpu-count {ncpu} '
           '--virtualbox-disk-size {diskmb} '
           '--virtualbox-memory {rammb} '
           '{name}').format(fs=args.storage_driver,
                            ncpu=args.numcpus,
                            diskmb=args.docker_machine_disk_size,
                            rammb=args.docker_machine_memory,
                            name=args.docker_machine_name)
    print cmd
    subprocess.check_call(cmd.split())


def _with_env(cmd):
    newcmd = ('eval $(docker-machine env %s) && ' % args.docker_machine_name) + cmd
    print '> ',newcmd
    return newcmd


def pull_image():
    img = image_name('moldesign_notebook')
    imgid = subprocess.check_output(_with_env('docker images -q %s' % img), shell=True).strip()
    if imgid == '':
        subprocess.check_call(_with_env('docker pull %s' % img), shell=True)


def compose_up(url, machine_name):
    print '\n_____________________________________\nStarting Buckyball and Jupyter ...'

    process = subprocess.Popen(_with_env('docker-compose -f %s up' ) % DOCKER_COMPOSE_PATH,
                               cwd=CONFIG_DIR, shell=True)
    atexit.register(down, 'stop')

    print 'Waiting for Jupyter server at %s ...' % url
    if not wait_net_service(url, JUPYTERPORT, timeout=20):
        raise IOError('Timeout waiting to connect to Jupyter @ %s:%s' % (url, JUPYTERPORT))
    return process


def compose_down(cmd):
    print 'Waiting for server to shut down cleanly ...'
    subprocess.check_output(_with_env('docker-compose -f %s %s' % (DOCKER_COMPOSE_PATH, cmd)),
                            cwd=CONFIG_DIR, shell=True)


URL_OPENERS = ['Open', 'xdg-open', 'sensible-browser', 'gnome-open', 'x-www-browser']


def open_jupyter(url):
    to_open = 'http://%s:%s' % (url, JUPYTERPORT)
    print '\n_____________________________________\nBuckyball is ready at\n\t\t %s' % to_open
    for exe in URL_OPENERS:
        if distutils.spawn.find_executable(exe) is not None:
            subprocess.check_call([exe, to_open])
            return
    print 'Point your browser to %s to get started.' % to_open  # fallback


is_mac = (check_path(['osascript'])['osascript'] is not None)


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


if __name__ == '__main__':
    args = make_parser().parse_args()
    main()
