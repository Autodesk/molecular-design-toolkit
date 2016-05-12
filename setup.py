import os
import sys
from os.path import relpath, join

import imp

import subprocess
from setuptools import find_packages, setup
from setuptools.command.install import install

assert sys.version_info[:2] == (2, 7), "Sorry, this package requires Python 2.7."

########################
__version__ = '0.3'
VERSION = __version__
ISRELEASED = False
########################
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

print 'printing works here'

class PostInstall(install):
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

setup(
    name='Molecular Design Toolkit',
    version='0.3',
    classifiers=CLASSIFIERS.split('\n'),
    packages=find_packages(),
    package_data={'moldesign': find_package_data()},
    install_requires=requirements,
    url='http://autodeskresearch.com',
    cmdclass={'install': PostInstall},
    license='Apache 2.0',
    author='Aaron Virshup',
    author_email='aaron.virshup [at] autodesk [dot] com',
    description='Dead-simple chemical simulation, visualization, and cloud computing in a notebook'
)
