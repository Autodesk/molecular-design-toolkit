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
This module stores definitions of common parameters for common techniques.

These are used to standardize our interfaces to other codes, and automatically generate interactive
notebook interfaces to configure various techniques.
"""
import operator as op

from . import units as u
from . import utils


def isin(a, b): return a in b


class WhenParam(object):
    def __init__(self, parameter, operator, checkval):
        self.operator = operator
        self.parameter = parameter
        self.checkval = checkval

    def __call__(self, paramset):
        """
        Args:
            paramset (dict):

        Returns:
            bool: True if the parameter is releveant, false otherwise
        """
        #TODO: anything relevant to an irrelevant parameter is also irrelevant
        return self.operator(paramset[self.parameter], self.checkval)


class Parameter(object):
    """ A generic parameter for a computational method

    Args:
        name (str): the arguments name (this is also its key in the method's ``params`` dictionary)
        short_description (str): A more readable description of about 100 characters
        type: The type of the param, including units if applicable.
            This may be a type (``int``, ``str``, etc.); if the quantity has physical units, you may
            also pass an example of this quantity (e.g., ``1.0 * units.angstrom``)
        default: the default value, or None if the user is required to set this parameter manually
        choices (list): A list of allowable values for the parameter
        help_url (str): URL for detailed help (not currently implemented)
        relevance (WhenParam): specifies when a given parameter will affect the dynamics

    Examples:
        >>> Parameter('timestep', 'Dynamics timestep', type=1.0*u.fs, default=2.0*u.fs)
        <Parameter "timestep", type: float, units: fs>
        >>> Parameter('functional', 'DFT XC functional', choices=['b3lyp', 'pbe0'],
        >>>           relevance=WhenParam('theory', op.eq, 'rks'))
        <Parameter "functional", type: str>
    """
    def __init__(self, name,
                 short_description=None,
                 type=None,
                 default=None,
                 choices=None,
                 help_url=None,
                 relevance=None):
        self.name = name
        self.displayname = utils.if_not_none(short_description, name)
        self.value = None
        self.default = default
        self.choices = utils.if_not_none(choices, [])
        self.type = type
        self.help_url = help_url
        if isinstance(type, u.MdtQuantity):
            type = type.units
        if isinstance(type, u.MdtUnit):
            self.type = float
            self.units = type
        else:
            self.units = None
        self.relevance = relevance

    def __str__(self):
        s = '%s "%s", type: %s' % (type(self).__name__, self.name, self.type.__name__)
        if self.units is not None:
            s += ', units: %s' % self.units
        return s

    def __repr__(self):
        try:
            return '<%s>' % self
        except (KeyError, AttributeError):
            return '<%s at %x - exception in __repr__>' % (type(self), id(self))



# TODO - make this ordered as well as dotted
def named_dict(l):
    return utils.DotDict((i.name, i) for i in l)

model_parameters = named_dict([
    Parameter('subsystem')
])

FORCEFIELDS = []
PERIODICITIES = [False, 'box']


mm_model_parameters = named_dict([
    Parameter('cutoff', 'Cutoff for nonbonded interactions',
              default=1.0*u.nm, type=u.nm),
    Parameter('nonbonded', 'Nonbonded interaction method', default='cutoff', type=str,
              choices=['cutoff', 'pme', 'ewald']),
    Parameter('implicit_solvent',
              'Implicit solvent method',
              type=str,
              choices=['gbsa', 'obc', 'pbsa', None]),
    Parameter('solute_dielectric', 'Solute dielectric constant',
              default=1.0, type=float),
    Parameter('solvent_dielectric', 'Solvent dielectric constant',
              default=78.5, type=float),
    Parameter('ewald_error', 'Ewald error tolerance', default=0.0005, type=float),
    Parameter('periodic', 'Periodicity', default=False, choices=PERIODICITIES)
])


QMTHEORIES = ['rhf', 'rks', 'mp2', 'casscf', 'casci', 'fci']
BASISSETS = ['3-21g', '4-31g', '6-31g', '6-31g*', '6-31g**',
             '6-311g', '6-311g*', '6-311g+', '6-311g*+',
             'sto-3g', 'sto-6g', 'minao', 'weigend',
             'dz' 'dzp', 'dtz', 'dqz',
             'aug-cc-pvdz', 'aug-cc-pvtz', 'aug-cc-pvqz']
FUNCTIONALS = ['b3lyp', 'blyp', 'pbe0', 'x3lyp', 'mpw3lyp5']
# This is a VERY limited set to start with; all hybrid functionals for now
# Need to think more about interface and what to offer by default
# PySCF xcs are at https://github.com/sunqm/pyscf/blob/master/dft/libxc.py for now

qm_model_parameters = named_dict([
    Parameter('theory', 'QM theory', choices=QMTHEORIES),
    Parameter('functional', 'DFT Functional', default='b3lyp',
              choices=FUNCTIONALS,  # TODO: allow separate x and c functionals
              relevance=WhenParam('theory', isin, 'dft rks ks uks'.split())),
    Parameter('active_electrons', 'Active electrons', type=int, default=2,
              relevance=WhenParam('theory', isin, ['casscf', 'mcscf', 'casci'])),
    Parameter('active_orbitals', 'Active orbitals', type=int, default=2,
              relevance=WhenParam('theory', isin, ['casscf', 'mcscf', 'casci'])),
    Parameter('state_average', 'States to average for SCF', type=int, default=1,
              relevance=WhenParam('theory', isin, ['casscf', 'mcscf'])),
    Parameter('basis', 'Basis set', choices=BASISSETS),
    Parameter('wfn_guess', 'Starting guess method', default='huckel',
              choices=['huckel', 'minao', 'stored']),
    Parameter('store_orb_guesses', 'Automatically use orbitals for next initial guess',
              default=True, type=bool),
    Parameter('multiplicity', 'Spin multiplicity', default=1, type=int),
    Parameter('symmetry', 'Symmetry detection',
              default=None, choices=[None, 'Auto', 'Loose']),
    Parameter('initial_guess', 'Wfn for initial guess',
              relevance=WhenParam('wfn_guess', op.eq, 'stored'))
              ])

integrator_parameters = named_dict([
    Parameter('timestep', 'Dynamics timestep', default=1.0*u.fs, type=u.default.time),
    Parameter('frame_interval', 'Time between frames', default=1.0*u.ps, type=u.fs)
])

md_parameters = named_dict([
    Parameter('remove_translation', 'Remove global translations', default=True,
              type=bool),
    Parameter('constrain_hbonds', 'Constrain covalent hydrogen bonds',
              default=True, type=bool),
    Parameter('constrain_water', 'Constrain water geometries',
              default=True, type=bool),
    Parameter('remove_rotation', 'Remove global rotations', default=False, type=bool),

])

constant_temp_parameters = named_dict([
    Parameter('temperature', 'Thermostat temperature', default=298 * u.kelvin,
              type=u.default.temperature)])

langevin_parameters = named_dict([
    Parameter('collision_rate', 'Thermal collision rate', default=1.0/u.ps, type=1/u.ps)
])

ground_state_properties = ['potential_energy',
                           'forces',
                           'dipole_moment',
                           'quadrupole_moment',
                           'octupole_moment',
                           'mulliken_charges',
                           'esp_charges',
                           'orbitals',
                           'orbital_energies',
                           'ci_vector',
                           'hessian',
                           'am1_bcc_charges']
"""If you're just calculating these, then just pass the
requested quantities as a list of keywords to the calculate method"""

excited_state_properties = ['state_energies',
                            'state_forces',
                            'state_ci_vector']
"""
When requesting these quantities, requests need to be passed to mol.calculate
as a dict with a list of states for each quantity, e.g.
>>> mol.calculate(requests={'state_energies':[1,2],'forces':[1,2]})
to get state_energies and forces for states 1 and 2.

Adiabatic states are indexed starting at 0, so state 0 is
the ground state, 1 is the first excited state, etc.
E.g.. state_energies[0] == potential_energy
"""

multistate_properties = ['transition_dipole',
                         'nacv',
                         'oscillator_strength']
"""
When requesting these quantities, requests need to be passed to mol.calculate
as a dict with a list of *pairs* of states for each quantity, e.g.
>>> mol.calculate(requests={'esp_charges':None, 'nacv':[(0,1),(0,2),(1,2)]})
"""
