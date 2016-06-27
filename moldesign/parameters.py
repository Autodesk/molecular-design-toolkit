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
This is where we store lists of keywords that should be consistent
across moldesign.

Consistency is not, and may never be, enforced.

TODO:
"""

from moldesign import units as u
from moldesign.utils import if_not_none, DotDict


class ForceField(object):
    """Generalized force field type (blank for now)"""
class BasisSet(object):
    """Generalized basis set type (blank for now"""
class ElectronicWfn(object):
    """Generalized orbital storage (blank for now)"""
class QMTheory(object):
    """Generalized QM type (blank for now)"""
class SymmetryGroup(object):
    """Generalized symmetry type (blank for now)"""


class Parameter(object):
    def __init__(self, name,
                 short_description=None,
                 type=None,
                 default=None,
                 choices=None,
                 select_multiple=False,
                 help_url=None):
        """
        A method's parameter
        :param default: the default value. If this does not match any spec in choices, it must be set.
        :param choices: A list of allowable values for the parameter
        :param types: a list of types for the parameter (does not check choices if parameter is one of these types)
        :param number: Number of values (>1 should be passed in list). '+' indicates arbitrary list size
        :return:
        """
        self.name = name
        self.value = None
        self.default = default
        self.choices = if_not_none(choices, [])
        self.type = type
        self.help_url = help_url
        if isinstance(type, u.MdtUnit):
            self.type = float
            self.units = type
        else:
            self.units = None
        self.select_multiple = select_multiple


def named_dict(l):
    return DotDict({i.name: i for i in l})

model_parameters = named_dict([
    Parameter('subsystem')
])

FORCEFIELDS = []
PERIODICITIES = [None, 'box']
QMTHEORIES = []
BASISSETS = []

mm_model_parameters = named_dict([
    Parameter('cutoff', 'Cutoff for nonbonded interactions', default=1.0*u.nm, type=u.nm),
    Parameter('nonbonded', 'Method to calculate nonbonded interactions', default='cutoff', type=str,
              choices=['cutoff', 'pme', 'ewald']),
    Parameter('implicit_solvent',
              'Implicit solvent method',
              type=str,
              choices=['gbsa', 'obc', 'pbsa', None]),
    Parameter('constrain_hbonds', 'Constrain covalent hydrogen bonds',
              default=True, type=bool),
    Parameter('constrain_water', 'Constrain water geometries',
              default=True, type=bool),
    Parameter('solute_dielectric', 'Solute dielectric constant',
              default=1.0, type=float),
    Parameter('solvent_dielectric', 'Solvent dielectric constant',
              default=78.5, type=float),
    Parameter('ewald_error', default=0.0005, type=float),
    Parameter('forcefield', type=ForceField,
              choices=FORCEFIELDS,
              select_multiple=True),
    Parameter('periodic', default=False, choices=PERIODICITIES)
])

qm_model_parameters = named_dict([
    Parameter('theory', choices=QMTHEORIES),
    Parameter('charge', default=None, type=int),  # not sure if this should be here or in mol
    Parameter('multiplicity', default=1, type=int),
    Parameter('basis_set', choices=BASISSETS),
    Parameter('symmetry', default=None, choices=[None, 'Auto', 'Loose']),
    Parameter('wfn_guess', default='huckel',
              choices=['huckel', 'guess']
              )])

integrator_parameters = named_dict([
    Parameter('timestep', default=1.0*u.fs, type=u.default.time),
    Parameter('remove_translation', default=True, type=bool),
    Parameter('remove_rotation', default=False, type=bool),
    Parameter('temperature', default=298 * u.kelvin, type=u.default.temperature),
    Parameter('collision_rate', default=1.0 / u.ps, type=[1.0 / u.ps]),
    Parameter('frame_interval', default=500, type=int)
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
