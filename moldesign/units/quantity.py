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

import operator
import copy
from os.path import join, abspath, dirname

import numpy as np
from pint import UnitRegistry, set_application_registry, DimensionalityError, UndefinedUnitError

from ..utils import ResizableArray


# Set up pint's unit definitions
ureg = UnitRegistry()
unit_def_file = join(abspath(dirname(__file__)), '..', '_static_data','pint_atomic_units.txt')
ureg.load_definitions(unit_def_file)
ureg.default_system = 'nano'
set_application_registry(ureg)


class MdtUnit(ureg.Unit):
    """
    Pickleable version of pint's Unit class.
    """
    def __reduce__(self):
        return _get_unit, (str(self),)

    @property
    def units(self):
        return self

    def convert(self, value):
        """ Returns quantity converted to these units

        Args:
            value (MdtQuantity or Numeric): value to convert

        Returns:
            MdtQuantity: converted value

        Raises:
            DimensionalityError: if the quantity does not have these units' dimensionality
        """
        if hasattr(value, 'to'):
            return value.to(self)
        elif self.dimensionless:
            return value * self
        else:
            raise DimensionalityError('Cannot convert "%s" to units of "%s"' % (value, self))

    def value_of(self, value):
        """ Returns numeric value of the quantity in these units

        Args:
            value (MdtQuantity or Numeric): value to convert

        Returns:
            Numeric: value in this object's units

        Raises:
            DimensionalityError: if the quantity does not have these units' dimensionality
        """
        v = self.convert(value)
        return v.magnitude


def _get_unit(unitname):
    """pickle helper for deserializing MdtUnit objects"""
    return getattr(ureg, unitname)


class MdtQuantity(ureg.Quantity):
    """
    This is a 'patched' version of pint's quantities that can be pickled (slightly hacky)
    and supports more numpy operations.
    Users should never need to instantiate this directly - instead, construct
    MDT quantities by multiplying numbers/arrays with the pre-defined units

    Examples:
        >>> 5.0 * units.femtoseconds
        >>> [1.0,2.0,3.0] * units.eV
    """
    # Patching some ufunc intercepts - these don't all necessarily work
    _Quantity__prod_units = ureg.Quantity._Quantity__prod_units.copy()
    _Quantity__prod_units['dot'] = 'mul'
    _Quantity__prod_units['cross'] = 'mul'

    _Quantity__copy_units = ureg.Quantity._Quantity__copy_units[:]
    _Quantity__copy_units.extend(('diagonal', 'append', '_broadcast_to'))
    _Quantity__handled = ureg.Quantity._Quantity__handled + ('diagonal', 'append', 'dot')

    # For pickling - prevent delegation to the built-in types' __getnewargs__ methods:
    def __getattr__(self, item):
        if item == '__getnewargs__':
            raise AttributeError('__getnewargs__ not accessible in this class')
        else:
            return super(MdtQuantity, self).__getattr__(item)

    def __reduce__(self):
        replacer = list(super(MdtQuantity, self).__reduce__())
        replacer[0] = MdtQuantity
        return tuple(replacer)

    def __deepcopy__(self, memo):
        result = copy.deepcopy(self.magnitude, memo) * self.get_units()
        memo[id(self)] = result
        return result

    def __hash__(self):
        m = self._magnitude
        if isinstance(m, np.ndarray) and m.shape == ():
            m = float(m)
        return hash((m, str(self.units)))

    def __setitem__(self, key, value):
        from . import array as quantityarray

        try:  # Speed up item assignment by overriding pint's implementation
            if self.units == value.units:
                self.magnitude[key] = value._magnitude
            else:
                self.magnitude[key] = value.value_in(self.units)
        except AttributeError:
            if isinstance(value, basestring):
                raise TypeError("Cannot assign units to a string ('%s')"%value)

            try:  # fallback to pint's implementation
                super().__setitem__(key, value)
            except (TypeError, ValueError):
                # one last ditch effort to create a more well-behaved object
                super().__setitem__(key, quantityarray(value))

    def __eq__(self, other):
        return self.compare(other, operator.eq)

    @property
    def shape(self):
        return self.magnitude.shape

    @shape.setter
    def shape(self, value):
        self.magnitude.shape = value

    def compare(self, other, op):
        """ Augments the :class:`pint._Quantity` method with the following features:
          - Comparisons to dimensionless 0 can proceed without unit checking
        """
        other = MdtQuantity(other)
        try:
            iszero = other.magnitude == 0.0 and other.dimensionless
        except ValueError:
            iszero = False

        if iszero:
            return op(self.magnitude, other.magnitude)
        else:
            return op(self.magnitude, other.value_in(self.units))

    def get_units(self):
        """
        Return the base unit system of an quantity
        """
        x = self
        while True:
            try: x = next(x.__iter__())
            except (AttributeError, TypeError): break
        try:
            y = 1.0 * x
            y._magnitude = 1.0
            return y
        except AttributeError:
            return 1.0

    def norm(self):
        """L2-norm of this object including units

        Returns:
            Scalar: L2-norm
        """
        units = self.get_units()
        return units * np.linalg.norm(self._magnitude)

    def normalized(self):
        """ Normalizes a vector or matrix

        Returns:
            np.ndarray: L2-normalized copy of this array (no units)
        """
        from ..mathutils import normalized
        return normalized(self.magnitude)

    def dot(self, other):
        """ Dot product that correctly multiplies units

        Returns:
            Array
        """
        if hasattr(other, 'get_units'):
            units = self.get_units() * other.get_units()
        else:
            units = self.get_units()
        return units * np.dot(self, other)

    def cross(self, other):
        """ Cross product that correctly multiplies units

        Returns:
            Array
        """

        if hasattr(other, 'get_units'):
            units = self.get_units() * other.get_units()
        else:
            units = self.get_units()
        return units * np.cross(self, other)

    def ldot(self, other):
        """ Left-multiplication version of dot that correctly multiplies units

        This is mathematically equivalent to ``other.dot(self)``, but preserves units even if ``other``
        is a plain numpy array

        Args:
            other (MdtQuantity or np.ndarray): quantity to take the dot product with

        Examples:
            >>> mat1 = np.ones((3,2))
            >>> vec1 = np.array([-3.0,2.0]) * u.angstrom
            >>> vec1.ldot(mat1)
            <Quantity([-1. -1. -1.], 'ang')>
            >>> # This won't work because "mat1", a numpy array, doesn't respect units
            >>> mat1.dot(vec1)
        """
        if hasattr(other, 'get_units'):
            units = self.get_units() * other.get_units()
        else:
            units = self.get_units()
        return units * np.dot(other, self)

    def __mod__(self, other):
        my_units = self.get_units()
        s_mag = self.magnitude
        o_mag = other.value_in(my_units)
        m = s_mag % o_mag
        return m * my_units

    # backwards-compatible name
    value_in = ureg.Quantity.m_as

    def defunits_value(self):
        return self.defunits().magnitude

    # defunits = ureg.Quantity.to_base_units  # replacing this with the new pint implementation
    def defunits(self):
        """Return this quantity in moldesign's default unit system (as specified in moldesign.units.default)"""
        from . import default
        return default.convert(self)

    # defunits_inplace = ureg.Quantity.ito_base_units  # replacing this with the new pint implementation
    def defunits_inplace(self):
        """Internally convert quantity to default units"""
        from . import default
        newunit = default.get_baseunit(self)
        return self.ito(newunit)

    def to_simtk(self):
        """ Return a SimTK quantity object
        """
        from moldesign.interfaces.openmm import pint2simtk
        return pint2simtk(self)

    def to_json(self):
        """ Convert to a simple JSON format

        Returns:
            dict: ``{value: <float>, units: <str>}``

        Examples:
            >>> from moldesign.units import angstrom
            >>> q = 1.0 * angstrom
            >>> q.to_json()
            {'units':'angstrom', value: 1.0}
        """
        mag = self.magnitude
        if isinstance(mag, np.ndarray):
            mag = mag.tolist()
        return {'value': mag,
                'units': str(self.units)}

    def make_resizable(self):
        self._magnitude = ResizableArray(self._magnitude)

    def append(self, item):
        mag = item.value_in(self.units)
        self._magnitude.append(mag)

    def extend(self, items):
        from . import array
        mags = array(items).value_in(self.units)
        self._magnitude.append(mags)

# monkeypatch pint's unit registry to return BuckyballQuantities
ureg.Quantity = MdtQuantity
ureg.Unit = MdtUnit

# These synonyms are here solely so that we can write descriptive docstrings
# TODO: use typing module to turn these into real abstract types, with dimensional parameterization

class Scalar(MdtQuantity):
    """ A scalar quantity (i.e., a single floating point number) with attached units
    """
    def __init__(self, *args):
        raise NotImplementedError('This is an abstract class - use MdtQuantity instead')


class Vector(MdtQuantity):
    """ A vector quantity (i.e., a list of floats) with attached units that behaves like a
    1-dimensional numpy array with units
    """
    def __init__(self, *args):
        raise NotImplementedError('This is an abstract class - use MdtQuantity instead')


class Array(MdtQuantity):
    """ A matrix quantity (i.e., a matrix of floats) with attached units that behaves like a
    2-dimensional numpy array with units
    """
    def __init__(self, *args):
        raise NotImplementedError('This is an abstract class - use MdtQuantity instead')


class Tensor(MdtQuantity):
    """ A multidimensional array of floats with attached units that behaves like a
    multidimensional numpy array with units
    """
    def __init__(self, *args):
        raise NotImplementedError('This is an abstract class - use MdtQuantity instead')

