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
    """
    Return the base unit system of an quantity
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
    try:
        y = 1.0 * x
        y._magnitude = 1.0
        return y
    except AttributeError:
        return 1.0 * ureg.dimensionless


def array(qlist, baseunit=None):
    """ Facilitates creating an array with units - like numpy.array, but it also checks
     units for all components of the array

     Args:
         qlist (List[MdtQuantity]): List-like object of quantity objects
         baseunit (MdtUnit) unit to standardize with

    Returns:
        MdtQuantity: array with standardized units
    """
    if hasattr(qlist, 'units') and hasattr(qlist, 'magnitude'):
        return MdtQuantity(qlist)

    if baseunit is None:
        baseunit = get_units(qlist)
        try:
            if baseunit == 1.0:
                return np.array(qlist)
        except DimensionalityError:
            pass

    try:
        newlist = [array(item, baseunit=baseunit).value_in(baseunit) for item in qlist]
        return baseunit * newlist
    except TypeError as exc:
        return qlist.to(baseunit)


#@utils.args_from(np.broadcast_to)
def broadcast_to(arr, *args, **kwargs):
    units = arr.units
    newarr = np.zeros(2) * units
    tmp = np.broadcast_to(arr, *args, **kwargs)
    newarr._magnitude = tmp
    return newarr
