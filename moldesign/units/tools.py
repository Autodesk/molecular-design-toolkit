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


def units_transfer(from_var, to_var, force=False):
    """ Give the "to_var" object the same units as "from_var"

    Args:
        from_var (MdtQuantity): use this quantities units
        to_var (MdtQuantity): apply units to this quantity
        force (bool):  Transfer the units even if from_var and to_var have incompatible units

    Returns:
        MdtQuantity: to_var with from_var's units
    """

    # If from_var is not a Quantity-like object, return a normal python scalar
    if not hasattr(from_var, '_units'):
        try:
            if to_var.dimensionless or force: return to_var._magnitude
            else: raise DimensionalityError(from_var, to_var)
        except AttributeError:
            return to_var

    # If to_var is a quantity-like object, just perform in-place conversion
    try:
        return to_var.to(from_var.units)
    except AttributeError:
        pass

    # If to_var has no units, return a Quantity object
    return to_var * from_var.units


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
    # TODO: replace with np.broadcast_to once numpy 1.10 is available in most package managers
    tmp = np.broadcast_to(arr, *args, **kwargs)
    newarr._magnitude = tmp
    return newarr


# TODO: remove once numpy 1.10 is available in most package managers
def _maybe_view_as_subclass(original_array, new_array):
    """ Temporary copy of a new NumPy function in 1.10, which isn't available in pkg managers yet
    """
    if type(original_array) is not type(new_array):
        # if input was an ndarray subclass and subclasses were OK,
        # then view the result as that subclass.
        new_array = new_array.view(type=type(original_array))
        # Since we have done something akin to a view from original_array, we
        # should let the subclass finalize (if it has it implemented, i.e., is
        # not None).
        if new_array.__array_finalize__:
            new_array.__array_finalize__(original_array)
    return new_array