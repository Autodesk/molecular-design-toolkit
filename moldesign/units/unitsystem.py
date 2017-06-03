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

from .constants import *


class UnitSystem(object):
    """ Class for standardizing units - specifies preferred units for length, mass, energy etc.

    In MDT, many methods will automatically convert output using the UnitSystem at
    ``moldesign.units.default``
    """
    def __init__(self, length, mass, time, energy,
                 temperature=kelvin,
                 force=None, momentum=None,
                 angle=radians,
                 charge=q_e):
        self.length = length
        self.mass = mass
        self.time = time
        self.energy = energy
        self.temperature = temperature
        self.force = force
        self.momentum = momentum
        self.angle = angle
        self.charge = charge

    def __getitem__(self, item):
        """ For convenience when using pint dimensionality descriptions.
        This aliases self['item'] = self['[item]'] = self.item,
        e.g. self['length'] = self['[length]'] = self.length
        """
        itemname = item.lstrip('[').rstrip(']')
        return getattr(self, itemname)

    @property
    def force(self):
        if self._force is None:
            return self.energy / self.length
        else:
            return self._force

    @force.setter
    def force(self, f):
        self._force = f

    @property
    def momentum(self):
        if self._momentum is None:
            return self.mass * self.length / self.time
        else:
            return self._momentum

    @momentum.setter
    def momentum(self, f):
        self._momentum = f

    def convert(self, quantity):
        """ Convert a quantity into this unit system.

        Args:
            quantity (MdtQuantity): quantity to convert
        """
        baseunit = self.get_baseunit(quantity)
        if isinstance(baseunit, int):
            assert baseunit == 1
            return quantity * ureg.dimensionless
        else:
            result = quantity.to(baseunit)
            return result

    def get_baseunit(self, quantity):
        """ Get units of a quantity, list or array

        Args:
            quantity (Any): any number or list-like object with units

        Raises:
            TypeError: if the passed object cannot have units (e.g., it's a string or ``None``)

        Returns:
            MdtUnit: units found in the passed object
        """
        try:
            dims = dict(quantity.dimensionality)
        except AttributeError:
            try:
                q = quantity[0]
            except (TypeError, StopIteration):
                raise TypeError('This type of object cannot have physical units')
            if isinstance(q, str):
                raise TypeError('This type of object cannot have physical units')
            try:
                return self.get_baseunit(q)
            except (IndexError, TypeError):  # Assume dimensionless
                return 1
        baseunit = 1

        # Factor out energy units (this doesn't work for things like energy/length)
        if '[length]' in dims and '[mass]' in dims and '[time]' in dims:
            while dims['[length]'] >= 1 and dims['[mass]'] >= 1 and dims['[time]'] <= -2:
                baseunit *= self['energy']
                dims['[length]'] -= 2
                dims['[mass]'] -= 1
                dims['[time]'] += 2

        if '[current]' in dims:
            dims.setdefault('[charge]', 0)
            dims.setdefault('[time]', 0)
            dims['[charge]'] += dims['[current]']
            dims['[time]'] -= dims['[current]']
            dims.pop('[current]')

        # Otherwise, just use the units
        for unit in dims:
            if dims[unit] == 0:
                continue
            try:
                baseunit *= self[unit]**dims[unit]
            except AttributeError:
                baseunit *= ureg[unit]**dims[unit]
        return baseunit


default = UnitSystem(length=angstrom, mass=amu, time=fs, energy=eV)
atomic_units = UnitSystem(length=a0, mass=m_e, time=t0, energy=hartree)
nano_si = UnitSystem(length=nm, mass=dalton, time=fs, energy=kjpermol)