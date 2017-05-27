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
"""
This module contains python "descriptors" (nothing to do with chemoinformatic "descriptors") that
help maintain the links between an atom's coordinates and its molecule's coordinates
"""

class ProtectedArray(object):
    """
    Descriptor for arrays that shouldn't be reassigned.
    Makes sure array attributes (specifically position and momentum) are modified in place

    Args:
        name (str): name of the instance attribute
    """
    def __init__(self, name):
        self.name = name

    def __get__(self, instance, cls=None):
        return getattr(instance, self.name)

    def __set__(self, instance, value):
        self.__get__(instance)[:] = value


class AtomArray(ProtectedArray):
    """
    Descriptor for atom coordinates are stored in the parent molecule.

    Makes sure that the arrays and their references are maintained during both
    reassignment and copying/pickling
    Args:
        atomname (str): name of the attribute in the atom instance
        parentname (str): name of the corresponding attribute in the molecule instance
    """
    def __init__(self, atomname, moleculename):
        self.name = atomname
        self.moleculename = moleculename

    def __get__(self, instance, cls=None):
        if instance.molecule is None:
            return getattr(instance, self.name)
        else:
            return getattr(instance.molecule, self.moleculename)[instance.index]


class AtomCoordinate(object):
    """ Descriptor for use with the ``Atom`` class.

    Gives access to 3D coordinates as ``atom.x,atom.y,atom.z`` instead of
    ``atom.position[0],atom.position[1],atom.position[2]``

    Args:
        quantity (str): name of the attribute that this accesses
        index (int): component of the attribute that this accesses
    """
    def __init__(self, attrname, index):
        self.attrname = attrname
        self.index = index

    def __get__(self, instance, cls):
        array = getattr(instance, self.attrname)
        return array[self.index]

    def __set__(self, instance, value):
        array = getattr(instance, self.attrname)
        array[self.index] = value