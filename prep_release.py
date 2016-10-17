#!/usr/bin/env python
""" This script builds and runs automated tests for a new moldesign release.
It's manual for now, will be replaced by Travis or Jenkins soon.

You shouldn't run this unless you know exactly what you're doing and why you're doing it
"""

import os
import sys
import subprocess
import atexit

version = sys.argv[1]


def run(cmd, display=True):
    print '\n>>> %s' % cmd
    if display:
        return subprocess.check_call(cmd, shell=True)
    else:
        return subprocess.check_output(cmd, shell=True)



DOCKER_REPO = 'autodesk/moldesign:'
DOCKER_RUN_SOCKET = 'docker run -v /var/run/docker.sock:/var/run/docker.sock'
COMPLETE_MDT = '%smoldesign_complete-%s' % (DOCKER_REPO, version)
MINIMAL_MDT = '%smoldesign_minimal-%s' % (DOCKER_REPO, version)
TEST_CMD = ('bash -c "pip install pytest pytest-xdist '
            "&& mkdir ~/.moldesign "
            "&& echo 'engine_type: docker' > ~/.moldesign/moldesign.yml "
            "&& cd /opt/molecular-design-toolkit/moldesign/_tests "
            '&& py.test -n 4"')


tags = set(run('git tag --list', display=False).splitlines())
rootdir = run('git rev-parse --show-toplevel', display=False).strip()

assert os.path.abspath(rootdir) == os.path.abspath(os.path.curdir), \
    "This command must be run at root of the repository directory"

# check that tag is valid
assert version not in tags, "Tag %s already exists!" % version
major, minor, patch = map(int, version.split('.'))


# Set the tag!
run('git tag %s' % version)
print 'Tag set: %s' % version


def untag():
    print 'Failed. Removing version tag %s' % version
    if not untag.success:
        run('git tag -d %s' % version)

untag.success = False

atexit.register(untag)

# Check that it propagated to the python package
import moldesign
assert os.path.abspath(os.path.join(moldesign.__path__[0],'..')) == os.path.abspath(rootdir)
assert moldesign.__version__ == version, 'Package has incorrect version: %s' % moldesign.__version__


# build docker images
run('cd docker_images '
    '&& ./docker-make.py --all --tag %s --repo autodesk/moldesign:' % version)


# Check dockerfile version tag
infolines = run(('docker run %s python -c '
                 '"import moldesign; print moldesign.__version__"') % COMPLETE_MDT,
                display=False).splitlines()
assert infolines[-1].strip() == version, \
    "moldesign in docker container reported wrong version: %s" % infolines[-1].strip()


# Check the remote image config version tag
infolines = run(('docker run %s python -c '
                 '"import moldesign; '
                 'print moldesign.compute.CONFIG_DEFAULTS.default_version_tag"') % COMPLETE_MDT,
                display=False).splitlines()
assert infolines[-1].strip() == version, \
    "moldesign.compute.config defaults to wrong image version: %s" % infolines[-1].strip()


# Run tests in "complete" image
run(' '.join([DOCKER_RUN_SOCKET, COMPLETE_MDT, TEST_CMD]))

# Run tests in "minimal" image
run(' '.join([DOCKER_RUN_SOCKET, MINIMAL_MDT, TEST_CMD]))

untag.success = True

# This is the irreversible part - so do it manually
print """
This LOOKS ready for release! Do the following to create version %s.
If you're ready, run these commands:
1. cd docker_images; ./docker-make.py --all --push --repo autodesk/moldesign: --tag %s
2. python setup.py register -r pypi
3. python setup.py sdist upload -r pypi
4. git push origin master --tags

Finally, mark it as release "v%s" in GitHub.

TO ABORT, RUN:
git tag -d %s
""" % (version, version, version, version)






