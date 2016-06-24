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
from moldesign import utils

from .quantity import *


def units_transfer(from_var, to_var, force=False):
    """
    Give the "to_var" object the same units as "from_var"
    :param from_var:
    :param to_var:
    :param force: Transfer the units even if from_var and to_var have incompatible units
    :return:
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
        try: x = x.__iter__().next()
        except (AttributeError, TypeError): break
    try:
        y = 1.0 * x
        y._magnitude = 1.0
        return y
    except AttributeError:
        return 1.0 * ureg.dimensionless


def array(qlist, baseunit=None):
    """ Facilitates creating an array with units - like numpy.array, but it also checks
     units for all components of the array

    :param qlist: List-like object of quantity objects
    :param baseunit: (optional) unit to standardize with
    :return: Quantity object
    """
    if baseunit is None:
        baseunit = get_units(qlist)
        if baseunit == 1.0: return np.array(qlist)

    try:
        newlist = [array(item, baseunit=baseunit).magnitude
                   for item in qlist]
        return baseunit * newlist
    except TypeError as exc:
        return qlist.to(baseunit)


@utils.args_from(np.broadcast_to)
def broadcast_to(arr, *args, **kwargs):
    units = arr.units
    newarr = np.zeros(2) * units
    tmp = np.broadcast_to(arr, *args, **kwargs)
    newarr._magnitude = tmp
    return newarr
