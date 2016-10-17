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

from setuptools import find_packages, setup
from setuptools.command.install import install

import versioneer

assert sys.version_info[:2] == (2, 7), "Sorry, this package requires Python 2.7."

PACKAGE_NAME = 'moldesign'

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


def find_package_data(pkgdir):
    """ Just include all files that won't be included as package modules.
    """
    files = []
    for root, dirnames, filenames in os.walk(pkgdir):
        not_a_package = '__init__.py' not in filenames
        for fn in filenames:
            basename, fext = os.path.splitext(fn)
            if not_a_package or (fext not in PYEXT) or ('static' in fn):
                files.append(relpath(join(root, fn), pkgdir))
    return files


class PostInstall(install):
    def run(self):
        install.run(self)
        self.prompt_intro()

    def prompt_intro(self):  # this doesn't actually display - print statements don't work?
        print 'Thank you for installing the Molecular Design Toolkit!!!'
        print 'For help, documentation, and any questions, visit us at '
        print '    http://moldesign.bionano.autodesk.com/'
        print '\nTo get started, please run:'
        print ' >>> python -m moldesign intro'

cmdclass = versioneer.get_cmdclass()
cmdclass['install'] = PostInstall

setup(
    name=PACKAGE_NAME,
    version=versioneer.get_version(),
    classifiers=CLASSIFIERS.splitlines(),
    packages=find_packages(),
    package_data={PACKAGE_NAME: find_package_data(PACKAGE_NAME)},
    install_requires=requirements,
    url='http://moldesign.bionano.autodesk.com',
    cmdclass=cmdclass,
    license='Apache 2.0',
    author='Aaron Virshup, BioNano Research at Atuodesk',
    author_email='moleculardesigntoolkit@autodesk.com',
    description='The Molecular Design Toolkit: Dead-simple chemical simulation, visualization, '
                'and cloud computing in a notebook'
)