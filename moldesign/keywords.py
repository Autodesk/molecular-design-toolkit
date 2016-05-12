"""
This is where we store lists of keywords that should be consistent
across moldesign.

Consistency is not, and may never be, enforced.

TODO:
"""
from moldesign.utils import if_not_none, DotDict
from moldesign import units as u

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
    def __init__(self,name,
                 default=None,
                 choices=None,
                 types=None,
                 number=1):
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
        self.choices = if_not_none(choices,[])
        self.types = if_not_none(types,[])


def named_dict(l):
    return DotDict({i.name: i for i in l})

model_parameters = named_dict([
    Parameter('subsystem')
])

mm_model_parameters = named_dict([
    Parameter('cutoff', default=1.0 * u.nm, types=u.nm, choices=[None]),
    Parameter('nonbonded', default='cutoff', choices=['cutoff', 'pme', 'ewald']),
    Parameter('implicit_solvent', choices=['gbsa', 'obc', 'pbsa', None]),
    Parameter('constrain_hbonds', default=True, types=[bool]),
    Parameter('constrain_water', types=[bool]),
    Parameter('solute_dielectric', default=1.0, types=[float]),
    Parameter('solvent_dielectric', default=78.5, types=[float]),
    Parameter('ewald_error', default=0.0005, types=[float, None]),
    Parameter('forcefield', types=[str, ForceField], number='+'),
    Parameter('periodic', default=False, types=[bool])])

qm_model_parameters = named_dict([
    Parameter('theory', types=[QMTheory]),
    Parameter('charge', default=None, types=[int]),
    Parameter('multiplicity', default=1, types=[int]),
    Parameter('basis_set', types=[BasisSet]),
    Parameter('symmetry', default=False, choices=[False, 'Auto', 'Loose'],
              types=[SymmetryGroup]),
    Parameter('wfn_guess', default='huckel',
              choices=['huckel', 'guess']
              )])

integrator_parameters = named_dict([
    Parameter('timestep', default=1.0 * u.fs, types=[u.fs]),
    Parameter('remove_translation', default=True, types=[bool]),
    Parameter('remove_rotation', default=False, types=[bool]),
    Parameter('temperature', default=298 * u.kelvin, types=[u.kelvin]),
    Parameter('collision_rate', default=1.0 / u.ps, types=[1.0 / u.ps]),
    Parameter('frame_interval', default=500, types=[int])
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
