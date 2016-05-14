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
"""
Set up physical constants and unit systems
"""
import copy
from os.path import join, abspath, dirname

import numpy as np
from pint import UnitRegistry, set_application_registry, DimensionalityError

# Set up pint's unit definitions
ureg = UnitRegistry()
unit_def_file = join(abspath(dirname(__file__)), 'static/pint_atomic_units.txt')
ureg.load_definitions(unit_def_file)
set_application_registry(ureg)


class MdtUnit(ureg.Unit):
    """
    Pickleable version of pint's Unit class.
    """
    def __reduce__(self):
        return get_unit, (str(self),)


def get_unit(unitname):
    return getattr(ureg, unitname)


class BuckyballQuantity(ureg.Quantity):
    """
    This is a 'patched' version of pint's quantities that can be pickled (slightly hacky)
    and supports more numpy operations.
    Users should never need to instantiate this directly - instead, construct
    BuckyBall quantities by multiplying numbers/arrays with the pre-defined units,
    e.g. 5.0 * units.femtoseconds or [1.0,2.0,3.0] * units.electronvolts
    """
    # Patching some ufunc intercepts - __prod_units does not seem to work, however
    _Quantity__prod_units = ureg.Quantity._Quantity__prod_units.copy()
    _Quantity__prod_units['dot'] = 'mul'
    _Quantity__prod_units['cross'] = 'mul'

    _Quantity__copy_units = ureg.Quantity._Quantity__copy_units[:]
    _Quantity__copy_units.extend(('diagonal', 'append'))
    _Quantity__handled = ureg.Quantity._Quantity__handled + ('diagonal', 'append', 'dot')

    # For pickling - prevent delegation to the built-in types' __getnewargs__ methods:
    def __getattr__(self, item):
        if item == '__getnewargs__':
            raise AttributeError('__getnewargs__ not accessible in this class')
        else:
            return super(BuckyballQuantity, self).__getattr__(item)

    def __reduce__(self):
        replacer = list(super(BuckyballQuantity, self).__reduce__())
        replacer[0] = BuckyballQuantity
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
    #        return super(BuckyballQuantity, self).__format__(fmt)

    def __setitem__(self, key, value):
        # Overrides pint's built-in version of this ... this is apparently way faster
        try:
            self.magnitude[key] = value.value_in(self.units)
        except AttributeError:
            if not hasattr(value, 'value_in'):  # deal with missing `value_in` method
                if self.dimensionless:  # case 1: this is OK if self is dimensionless
                    self.magnitude[key] = value
                else:  # case 2: User tried to pass a number without units
                    raise DimensionalityError('%s cannot be assigned to array with dimensions %s' %
                                              (value, self.units))
            else:  # case 3: attribute error is unrelated to this
                raise



    def __eq__(self, other):
        # override pint's implementation, which appears to give inconsistent results?
        if hasattr(other, 'value_in'):
            return self.magnitude == other.value_in(self.units)
        elif self.dimensionless:  # they're both dimensionless
            return self.magnitude == other
        else:  # only 'other' is dimensionless
            return False

    def get_units(self):
        """
        Return the base unit system of an quantity
        """
        x = self
        while True:
            try: x = x.__iter__().next()
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

    def dot(self, other):
        """Compute norm but respect units
        Returns:
            BuckyballQuantity
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

    def value_in(self, units):
        val = self.to(units)
        return val._magnitude

    def defunits_value(self):
        return self.defunits()._magnitude

    # defunits = ureg.Quantity.to_base_units  # replacing this with the new pint implementation
    def defunits(self):
        """Return this quantity in moldesign's default unit system (as specified in moldesign.units.default)"""
        import moldesign.units as u
        return u.default.convert(self)

    # defunits_inplace = ureg.Quantity.ito_base_units  # replacing this with the new pint implementation
    def defunits_inplace(self):
        """Internally convert quantity to default units"""
        import moldesign.units as u
        newunit = u.default.get_baseunit(self)
        return self.ito(newunit)

    def to_simtk(self):
        from moldesign.interfaces.openmm import pint2simtk
        return pint2simtk(self)

# monkeypatch pint's unit registry to return BuckyballQuantities
ureg.Quantity = BuckyballQuantity
ureg.Unit = MdtUnit

# These synonyms are here solely so that we can write descriptive docstrings
Scalar = Vector = Array = Tensor = BuckyballQuantity

# Constants
unity = ureg.angstrom / ureg.angstrom
imi = 1.0j
pi = np.pi
sqrtpi = np.sqrt(np.pi)
sqrt2 = np.sqrt(2.0)
epsilon_0 = ureg.epsilon_0
c = ureg.speed_of_light
alpha = ureg.fine_structure_constant
hbar = ureg.hbar
boltz = boltzmann_constant = k_b = ureg.boltzmann_constant
avogadro = (1.0 * ureg.mole * ureg.avogadro_number).to(unity).magnitude

# atomic units
hartree = ureg.hartree
a0 = bohr = ureg.bohr
atomic_time = t0 = ureg.t0
electron_mass = m_e = ureg.electron_mass
electron_charge = q_e = ureg.elementary_charge

# useful units
fs = femtoseconds = ureg.fs
ps = picoseconds = ureg.ps
eV = electronvolts = ureg.eV
kcalpermol = ureg.kcalpermol
gpermol = ureg.gpermol
kjpermol = ureg.kjpermol
radians = rad = ureg.rad
degrees = deg = ureg.degrees
amu = da = dalton = ureg.amu
kelvin = ureg.kelvin
nm = ureg.nanometers
ang = angstrom = ureg.ang

# sets default unit systems
def_length = angstrom
def_time = fs
def_vel = angstrom / fs
def_mass = amu
def_momentum = def_mass * def_vel
def_force = def_momentum / def_time
def_energy = eV


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
        return 1.0


def to_units_array(qlist, baseunit=None):
    """
    Facilitates creating a numpy array by standardizing units across the entire object
    :param qlist: List-like object of quantity objects
    :param baseunit: (optional) unit to standardize with
    :return: Quantity object
    """
    if baseunit is None:
        baseunit = get_units(qlist)
        if baseunit == 1.0: return np.array(qlist)

    try:
        newlist = [to_units_array(item,baseunit=baseunit).magnitude
                   for item in qlist]
        return baseunit * newlist
    except TypeError as exc:
        return qlist.to(baseunit)


class UnitSystem(object):
    """
    Class for standardizing units
    Many methods will look for a unit system at moldesign.units.default
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
        """
        Convert a quantity into this unit system
        @param quantity: moldesign.external.pint.Quantity
        """
        baseunit = self.get_baseunit(quantity)
        result = quantity.to(baseunit)
        return result

    def get_baseunit(self, quantity):
        # TODO: this needs to deal with angles
        # if quantity.dimensionless: return 1.0 # don't call this - it's super-slow (Pint 0.6.1)

        try:
            dims = dict(quantity.dimensionality)
        except AttributeError:
            try: return self.get_baseunit(quantity[0])
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
            if dims[unit] == 0: continue
            try:
                baseunit *= self[unit]**dims[unit]
            except AttributeError:
                baseunit *= ureg[unit]**dims[unit]
        return baseunit


default = UnitSystem(length=angstrom, mass=amu, time=fs, energy=eV)
atomic_units = UnitSystem(length=a0, mass=m_e, time=t0, energy=hartree)
nano_si = UnitSystem(length=nm, mass=dalton, time=fs, energy=kjpermol)