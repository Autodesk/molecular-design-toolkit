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
from os.path import relpath, join

import imp

import subprocess
from setuptools import find_packages, setup
from setuptools.command.install import install

import versioneer
cmdclass = versioneer.get_cmdclass()

assert sys.version_info[:2] == (2, 7), "Sorry, this package requires Python 2.7."


CLASSIFIERS = """\
Development Status :: 4 - Beta
Intended Audience :: Science/Research
Intended Audience :: Developers
Intended Audience :: Education
License :: OSI Approved :: Apache Software License
Programming Language :: Python :: 2.7
Programming Language :: Python :: 2 :: Only
Topic :: Scientific/Engineering :: Chemistry
Topic :: Scientific/Engineering :: Physics
Topic :: Scientific/Engineering :: Visualization
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

HOME = os.environ['HOME']
CONFIG_DIR = os.path.join(HOME, '.moldesign')
PYEXT = set('.py .pyc .pyo'.split())

with open('requirements.txt', 'r') as reqfile:
    requirements = [x.strip() for x in reqfile if x.strip()]

def find_package_data():
    files = []
    for root, dirnames, filenames in os.walk('moldesign'):
        for fn in filenames:
            basename, fext = os.path.splitext(fn)
            if fext not in PYEXT or ('static' in fn):
                files.append(relpath(join(root, fn), 'moldesign'))
    return files


class PostInstall(cmdclass['install']):
    def run(self):
        install.run(self)
        #self.protect_db()  # causes pip install to crash
        print 'hi'
        self.prompt_intro()

    def protect_db(self):
        # Prevent residue dictionary from getting corrupted
        print 'heyhey'
        modpath = imp.find_module('moldesign')[1]
        dbpath = os.path.join(modpath, 'static/residue_bonds')
        subprocess.check_call('chmod a-w {0}.dir {0}.dat'.format(dbpath))

    def prompt_intro(self):  # this doesn't actually display - print statements don't work?
        print 'Thank you for installing the Molecular Design Toolkit!!!'
        print 'For help, documentation, and any questions, visit us at '
        print '    http://bionanoresearch.com/moldesign'
        print '\nTo get started, please run:'
        print ' >>> python -m moldesign intro'

cmdclass['install'] = PostInstall

setup(
    name='Molecular Design Toolkit',
    version=versioneer.get_version(),
    classifiers=CLASSIFIERS.split('\n'),
    packages=find_packages(),
    package_data={'moldesign': find_package_data()},
    install_requires=requirements,
    url='http://autodeskresearch.com',
    cmdclass=cmdclass,
    license='Apache 2.0',
    author='Aaron Virshup',
    author_email='aaron.virshup [at] autodesk [dot] com',
    description='Dead-simple chemical simulation, visualization, and cloud computing in a notebook'
)