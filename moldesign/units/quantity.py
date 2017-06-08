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
import operator
import copy
from os.path import join, abspath, dirname
import numbers

import numpy as np
from pint import UnitRegistry, set_application_registry, DimensionalityError

from ..utils import ResizableArray

# Set up pint's unit definitions
ureg = UnitRegistry()
unit_def_file = join(abspath(dirname(__file__)), '../_static_data/pint_atomic_units.txt')
ureg.load_definitions(unit_def_file)
set_application_registry(ureg)


class MdtUnit(ureg.Unit):
    """
    Pickleable version of pint's Unit class.
    """
    def __reduce__(self):
        return _get_unit, (str(self),)


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

    # This doesn't deal with length specs correctly (pint's doesn't either though)
    #def __format__(self, fmt):
    #    fmtstring = '{m:%s} {u}' % fmt
    #    try:
    #        return fmtstring.format(m=self.magnitude,
    #                                u=self.units)
    #    except:
    #        return super().__format__(fmt)

    def __setitem__(self, key, value):
        # Overrides pint's built-in version of this ... this is apparently way faster
        try:
            self.magnitude[key] = value.value_in(self.units)
        except AttributeError:
            if not hasattr(value, 'value_in'):  # deal with missing `value_in` method
                if self.dimensionless:  # case 1: this is OK if self is dimensionless
                    self.magnitude[key] = value
                elif not isinstance(value, numbers.Number):  # case 2: this is not a number
                    raise TypeError('"%s" is not a valid numeric value' % value)
                else:  # case 3: wrong units
                    raise DimensionalityError(self.units, ureg.dimensionless)
            else:  # case 3: attribute error is unrelated to this
                raise

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
        """Compute norm but respect units"""
        units = self.get_units()
        return units * np.linalg.norm(self._magnitude)

    def normalized(self):
        return self/self.norm()

    def dot(self, other):
        """ Dot product that correctly multiplies units

        Returns:
            MdtQuantity
        """
        if hasattr(other, 'get_units'):
            units = self.get_units() * other.get_units()
        else:
            units = self.get_units()
        return units * np.dot(self, other)

    def cross(self, other):
        if hasattr(other, 'get_units'):
            units = self.get_units() * other.get_units()
        else:
            units = self.get_units()
        return units * np.cross(self, other)

    def ldot(self, other):
        """
        Left-multiplication version of dot.
        Use this to preserve units (built-in numpy versions don't)
        getting hackier ...
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
# TODO: use typing module to turn these into real abstract types

class Scalar(MdtQuantity):
    """ A scalar quantity (i.e., a single floating point number) with attached units
    """
    def __init__(self, *args):
        raise NotImplementedError('This is an abstract class - use MdtQuantity instead')


class Vector(MdtQuantity):
    """ A vector quantity (i.e., a list of floats) with attached units, which behaves like a
    1-dimensional numpy array
    """
    def __init__(self, *args):
        raise NotImplementedError('This is an abstract class - use MdtQuantity instead')


class Array(MdtQuantity):
    """ A matrix quantity (i.e., a matrix of floats) with attached units, which behaves like a
    2-dimensional numpy array
    """
    def __init__(self, *args):
        raise NotImplementedError('This is an abstract class - use MdtQuantity instead')


class Tensor(MdtQuantity):
    """ A vector quantity (i.e., a list of floats) with attached units, which behaves like a
    multidimensional numpy array
    """
    def __init__(self, *args):
        raise NotImplementedError('This is an abstract class - use MdtQuantity instead')

