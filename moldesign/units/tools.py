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

from .. import utils

from .quantity import *



def unitsum(iterable):
    """
    Faster method to compute sums of iterables if they're all in the right units

    Args:
        iterable (Iter[MdtQuantity]): iterable to sum over

    Returns:
        MdtQuantity: the sum
    """
    g0 = next(iterable).copy()
    for item in iterable:
        if item.units == g0.units:
            g0._magnitude += item._magnitude
        else:
            g0 += item
    return g0


def dot(a1, a2):
    """ Dot product that respects units

    Args:
        a1 (MdtQuantity or np.ndarray): First term in dot product
        a2 (MdtQuantity or np.ndarray): Second term in dot product

    Returns:
        MdtQuantity or np.ndarray: dot product (MdtQuantity if either input has units, ndarray else)
    """
    if isinstance(a2, MdtQuantity):
        return a2.ldot(a1)
    else:  # this will work whether or not a1 has units
        return a1.dot(a2)


@utils.args_from(np.linspace)
def linspace(start, stop, **kwargs):
    u1 = getattr(start, 'units', ureg.dimensionless)
    u2 = getattr(stop, 'units', ureg.dimensionless)
    if u1 == u2 == ureg.dimensionless:
        return np.linspace(start, stop, **kwargs)
    else:
        q1mag = start.magnitude
        q2mag = stop.value_in(start.units)
        return np.linspace(q1mag, q2mag, **kwargs) * start.units


def arrays_almost_equal(a1, a2):
    """ Return true if arrays are almost equal up to numerical noise

    Note:
        This is assumes that absolute differences less than 1e-12 are insignificant. It is
        therefore more likely to return "True" for very small numbers and
        "False" for very big numbers. Caveat emptor.

    Args:
        a1 (MdtQuantity or np.ndarray): first array
        a2 (MdtQuantity or np.ndarray): second array

    Returns:
        bool: True if arrays differ by no more than numerical noise in any element

    Raises:
        DimensionalityError: if the arrays have incompatible units
    """

    a1units = False
    if isinstance(a1, MdtQuantity):
        if a1.dimensionless:
            a1mag = a1.value_in(ureg.dimensionless)
        else:
            a1units = True
            a1mag = a1.magnitude

    else:
        a1mag = a1

    if isinstance(a2, MdtQuantity):
        if a2.dimensionless:
            if a1units:
                raise DimensionalityError(a1.units, ureg.dimensionless,
                                          "Cannot compare objects")
            else:
                a2mag = a2.value_in(ureg.dimensionless)
        elif not a1units:
            raise DimensionalityError(ureg.dimensionless, a2.units,
                                      "Cannot compare objects")
        else:
            a2mag = a2.value_in(a1.units)
    else:
        if a1units:
            raise DimensionalityError(a1.units, ureg.dimensionless,
                                      "Cannot compare objects")
        else:
            a2mag = a2

    return np.allclose(a1mag, a2mag, atol=1e-12)


def from_json(j):
    """
    Convert a JSON description to a quantity.
    This is the inverse of :meth:`moldesign.units.quantity.MDTQuantity.to_json`

    Args:
        j (dict): ``{value: <float>, units: <str>}``

    Returns:
        moldesign.units.quantity.MDTQuantity
    """
    return j['value'] * ureg(j['units'])


def get_units(q):
    """ Return the base unit system of an quantity or arbitrarily-nested iterables of quantities

    Note: This routine will dive on the first element of iterables until a quantity with units
      until the units can be determined. It will not check the remaining elements of the iterable
      for consistency

    Examples:
        >>> from moldesign import units
        >>> units.get_units(1.0 * units.angstrom)
        <Unit('angstrom')>
        >>> units.get_units(np.array([1.0, 2, 3.0]))
        <Unit('dimensionless')>
        >>> # We dive on the first element of each iterable until we can determine a unit system:
        >>> units.get_units([[1.0 * u.dalton, 3.0 * u.eV], ['a'], 'gorilla'])
        <Unit('amu')>

    Args:
        q (MdtQuantity or numeric): quantity to test

    Returns:
        MdtUnit: the quantity's units
    """
    x = q
    while True:
        try:
            x = next(x.__iter__())
        except (AttributeError, TypeError):
            break
        else:
            if isinstance(x, str):
                raise TypeError('Found string data while trying to determine units')

    q = MdtQuantity(x)
    if q.dimensionless:
        return ureg.dimensionless
    else:
        return q.units


def array(qlist, defunits=False, _baseunit=None):
    """ Facilitates creating an array with units - like numpy.array, but it also checks
     units for all components of the array.

     Note:
         Unlike numpy.array, these arrays must have numeric type - this routine will
         raise a ValueError if a non-square array is passed.

     Args:
         qlist (List[MdtQuantity]): List-like object of quantity objects
         defunits (bool): if True, convert the array to the default units

    Returns:
        MdtQuantity: ndarray-like object with standardized units

    Raises:
        DimensionalityError: if the array has inconsistent units
        ValueError: if the array could not be converted to a square numpy array
    """
    from . import default

    if hasattr(qlist, 'units') and hasattr(qlist, 'magnitude'):
        return MdtQuantity(qlist)

    if _baseunit is None:
        _baseunit = get_units(qlist)
        if _baseunit.dimensionless:
            return _make_nparray(qlist)
        if defunits:
            _baseunit = default.get_default(_baseunit)

    if hasattr(qlist, 'to'):  # if already a quantity, just convert and return
        return qlist.to(_baseunit)

    try:  # try to create a quantity
        return _baseunit * [array(item, _baseunit=_baseunit).value_in(_baseunit) for item in qlist]
    except TypeError:  # if here, one or more objects cannot be converted to quantities
        raise DimensionalityError(_baseunit, ureg.dimensionless,
                                  extra_msg='Object "%s" does not have units' % qlist)


def _make_nparray(q):
    """ Turns a list of dimensionless numbers into a numpy array. Does not permit object arrays
    """
    if hasattr(q, 'units'):
        return q.value_in(ureg.dimensionless)
    try:
        arr = np.array([_make_nparray(x) for x in q])
        if arr.dtype == 'O':
            raise ValueError("Could not create numpy array of numeric data - is your input square?")
        else:
            return arr
    except TypeError:
        return q


#@utils.args_from(np.broadcast_to)
def broadcast_to(arr, *args, **kwargs):
    units = arr.units
    newarr = np.zeros(2) * units
    tmp = np.broadcast_to(arr, *args, **kwargs)
    newarr._magnitude = tmp
    return newarr
