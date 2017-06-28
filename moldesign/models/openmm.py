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
from future.utils import native_str

import moldesign.molecules
from moldesign import compute
from ..molecules import Trajectory, MolecularProperties
from ..utils import exports
from ..interfaces import openmm as opm
from .base import MMBase
from .. import parameters

openmm_platform_selector = parameters.Parameter(
        'compute_platform', 'OpenMM computing platform',
        type=str, default='cpu', choices=['opencl', 'cuda', 'cpu', 'reference', 'auto'],
        help_url='http://docs.openmm.org/7.1.0/userguide/library.html#platforms')

numcpus = parameters.num_cpus.copy()
numcpus.relevance = parameters.WhenParam('compute_platform', parameters.op.eq, 'cpu')
numcpus.help = ('Sets number of threads for OpenMM CPU platform. '
                'If 0, uses OpenMM defaults (which can be controlled '
                'via the OPENMM_NUM_THREADS environment variable). ')
numcpus.help_url = \
    'http://docs.openmm.org/7.1.0/userguide/library.html#platform-specific-properties'


@exports
class OpenMMPotential(MMBase, opm.OpenMMPickleMixin):
    """Creates an OpenMM "context" to drive energy calculations.
    Note that, while a dummy integrator is assigned, a different context will
    be created for any MD calculations.

    Attributes:
       sim (simtk.openmm.app.Simulation): OpenMM simulation object (once created)
    """
    # NEWFEATURE: need to set/get platform (and properties, e.g. number of threads)
    DEFAULT_PROPERTIES = ['potential_energy', 'forces']
    PARAMETERS = MMBase.PARAMETERS + [openmm_platform_selector, numcpus]
    _CALLS_MDT_IN_DOCKER = opm.force_remote

    _openmm_compatible = True

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._reset()

    def _reset(self):
        self.sim = None
        self.mm_system = None
        self.mm_integrator = None
        self._prepped_integrator = 'uninitialized'
        self._constraints_set = False
        self._required_tolerance = None

    def get_openmm_simulation(self):
        if opm.force_remote:
            raise ImportError("Can't create an OpenMM object on this machine - OpenMM not "
                              "installed")
        else:
            if not self._prepped: self.prep()
            return self.sim

    @compute.runsremotely(enable=opm.force_remote, is_imethod=True)
    def calculate(self, requests=None):
        """
        Drive a calculation and, when finished, update the parent molecule with the results.
        TODO: this update is SYNCHRONOUS, unlike other `calculate` methods that run remotely.
        TODO: Probably need to update DummyJob (call it localjob, syncjob?) to handle this
        :param requests: list of quantities to calculate
        :return: PythonJob-like object
        """
        if requests is None:
            requests = self.DEFAULT_PROPERTIES
        self.prep()
        self._set_openmm_state()
        state = self.sim.context.getState(getForces=True, getEnergy=True)
        props = MolecularProperties(self.mol,
                                    potential_energy=opm.simtk2pint(state.getPotentialEnergy()),
                                    forces=opm.simtk2pint(state.getForces(), flat=False))
        return props

    def prep(self):
        """ Construct the OpenMM simulation objects

        Note:
            An OpenMM simulation object consists of both the system AND the integrator. This routine
            therefore constructs both. If self.mol does not use an OpenMM integrator, we create
            an OpenMM simulation with a "Dummy" integrator that doesn't ever get used.
        """
        if opm.force_remote:
            return True

        from simtk.openmm import app

        if getattr(self.mol.integrator, '_openmm_compatible', False):
            _prepped = (self._prepped and self.mol.integrator._prepped and
                        self.mol.integrator is self._prepped_integrator)
            setup_integrator = True
        else:
            _prepped = self._prepped
            setup_integrator = False

        if _prepped:
            return

        self._reset()
        system_params = self._get_system_params()
        self.mm_system = self.mol.ff.parmed_obj.createSystem(**system_params)

        if setup_integrator:
            try:
                self._set_constraints()
            except moldesign.NotSupportedError as exc:
                print("Warning: dynamics not supported: %s" % exc.args[0])
            self.mm_integrator = self.mol.integrator.get_openmm_integrator()
        else:
            self.mm_integrator = self._make_dummy_integrator()

        platform, platform_properties = self._get_platform()

        if self._required_tolerance:
            self.mm_integrator.setConstraintTolerance(float(self._required_tolerance))
        self.sim = app.Simulation(self.mol.ff.parmed_obj.topology,
                                  self.mm_system,
                                  self.mm_integrator,
                                  platform=platform,
                                  platformProperties=platform_properties)

        if setup_integrator:
            self.mol.integrator.energy_model = self
            self.mol.integrator.sim = self.sim
            self.mol.integrator._prepped = True

        self._prepped = True
        self._prepped_integrator = self.mol.integrator
        print('Created OpenMM kernel (Platform: %s)' % self.sim.context.getPlatform().getName())

    def minimize(self, **kwargs):
        if self.constraints_supported():
            traj = self._minimize(**kwargs)
            if opm.force_remote or (not kwargs.get('wait', False)):
                self._sync_remote(traj.mol)
                traj.mol = self.mol
            return traj
        else:
            return super().minimize(**kwargs)

    def _sync_remote(self, mol):
        # TODO: this is a hack to update the molecule object after a minimization
        #        We need a better pattern for this, ideally one that doesn't
        #        require an explicit wrapper like this - we shouldn't have to copy
        #        the properties over manually
        self.mol.positions = mol.positions
        self.mol.momenta = mol.momenta
        self.mol.properties = mol.properties
        self.mol.time = mol.time

    @compute.runsremotely(enable=opm.force_remote, is_imethod=True)
    def _minimize(self, nsteps=500,
                 force_tolerance=None,
                 frame_interval=None):
        """ Run an OpenMM minimization.
        Note:
            We're not able to support `frame_interval` for this method;
                see https://github.com/pandegroup/openmm/issues/1155
        Args:
           nsteps(int): maximum number of steps
           force_tolerance (moldesign.units.MdtQuantity): RMSD tolerance for convergence
               [energy/length]
        """

        # NEWFEATURE: write/find an openmm "integrator" to do this minimization.
        # openmm doesn't work with small step numbers, and doesn't support
        # callbacks during minimization, so frame_interval is disabled.
        if frame_interval is not None:
            raise ValueError('frame_interval not supported by OpenMM minimizer.'
                             ' Use a method from moldesign.minimizers instead')

        self.prep()
        trajectory = Trajectory(self.mol)
        self._set_openmm_state()
        trajectory.new_frame(annotation='initial structure, energy=%s' %
                                        self.mol.calc_potential_energy().defunits())

        self.sim.minimizeEnergy(maxIterations=nsteps, tolerance=0.01)  # tolerance is in kj/mole ...
        self._sync_to_openmm()
        new_energy = self.mol.calc_potential_energy()
        trajectory.new_frame(annotation='minimization result, energy=%s' % new_energy.defunits())

        return trajectory

    def constraints_supported(self):
        """ Check whether this molecule's constraints can be enforced in OpenMM

        This sort of overlaps with _set_constraints, but doesn't have dependencies on OpenMM
        """
        for constraint in self.mol.constraints:
            if constraint.desc not in ('position', 'distance', 'hbonds'):
                return False
        else:
            return True

    def _set_constraints(self):
        if self._constraints_set:
            return
        system = self.mm_system
        fixed_atoms = set()

        # openmm uses a global tolerance, calculated as ``constraint_violation/constraint_dist``
        # (since only distance constraints are supported). Here we calculate the necessary value
        required_tolerance = None

        # Constrain atom positions
        for constraint in self.mol.constraints:
            if constraint.desc == 'position':
                fixed_atoms.add(constraint.atom)
                self.mol.assert_atom(constraint.atom)
                system.setParticleMass(constraint.atom.index, 0.0)

            # Constrain distances between atom pairs
            elif constraint.desc == 'distance':
                self.mol.assert_atom(constraint.a1)
                self.mol.assert_atom(constraint.a2)
                system.addConstraint(constraint.a1.index,
                                     constraint.a2.index,
                                     opm.pint2simtk(constraint.value))

                if required_tolerance is None:
                    required_tolerance = 1e-5
                required_tolerance = min(required_tolerance, constraint.tolerance/constraint.value)

            elif constraint.desc == 'hbonds':
                continue  # already dealt with at system creation time

            else:
                raise moldesign.NotSupportedError("OpenMM does not support '%s' constraints" %
                                                  constraint.desc)

        # Workaround for OpenMM issue: can't have an atom that's both
        # fixed *and* has a distance constraint. If both atoms in the distance constraint are
        # also fixed, then we can just remove the constraint
        if len(fixed_atoms) > 0:
            num_constraints = system.getNumConstraints()
            ic = 0
            while ic < num_constraints:
                i, j, dist = system.getConstraintParameters(ic)
                ai = self.mol.atoms[i]
                aj = self.mol.atoms[j]
                if (ai in fixed_atoms) and (aj in fixed_atoms):
                    system.removeConstraint(ic)
                    num_constraints -= 1
                elif (ai in fixed_atoms) or (aj in fixed_atoms):  # only one is fixed
                    raise ValueError('In OpenMM, fixed atoms cannot be part of a constrained '
                                     'bond (%s)' % moldesign.molecules.bonds.Bond(ai, aj))
                else:
                    ic += 1

        self._constraints_set = True
        self._required_tolerance = required_tolerance

    @staticmethod
    def _make_dummy_integrator():
        from simtk import unit as stku
        from simtk import openmm
        return openmm.VerletIntegrator(2.0 * stku.femtoseconds)

    def _set_openmm_state(self):  # TODO: periodic state
        self.sim.context.setPositions(opm.pint2simtk(self.mol.positions))
        self.sim.context.setVelocities(opm.pint2simtk(self.mol.velocities))
        self.sim.context.setTime(opm.pint2simtk(self.mol.time))

    def _sync_to_openmm(self, positions=True, momenta=True, time=True):
        """
        Syncs the moldesign molecule's positions, momenta, and time to the simulation's
        """
        state = None
        get_positions = (positions is True)
        get_velocities = (momenta is True)
        if get_positions or get_velocities:
            state = self.sim.context.getState(getPositions=get_positions,
                                              getVelocities=get_velocities)
            if get_positions:
                positions = state.getPositions()
            if get_velocities:
                velocities = state.getVelocities()

            for iatom, atom in enumerate(self.mol.atoms):
                if positions != False: #this is weird because of numpy comparisons
                    atom.position = opm.simtk2pint(positions[iatom])
                if velocities != False:
                    atom.momentum = opm.simtk2pint(velocities[iatom]) * atom.mass

        if time is True:
            if state is None:
                state = self.sim.context.getState()
            time = opm.simtk2pint(state.getTime())
        elif time:
            time = opm.simtk2pint(time)

        if time:
            self.mol.time = time

    def _get_system_params(self):
        """ Translates the spec from MMBase into system parameter keywords for createSystem
        """
        # need cmm motion
        from simtk.openmm import app
        nonbonded_names = {'nocutoff': app.NoCutoff,
                           'ewald': app.Ewald,
                           'pme': app.PME,
                           'cutoff': (app.CutoffPeriodicif
                                      if self.params.periodic
                                      else app.CutoffNonPeriodic)}
        implicit_solvent_names = {'obc': app.OBC2,
                                  'obc1': app.OBC1,
                                  None: None}

        system_params = dict(nonbondedMethod=nonbonded_names[self.params.nonbonded],
                             nonbondedCutoff=opm.pint2simtk(self.params.cutoff),
                             implicitSolvent=implicit_solvent_names[self.params.implicit_solvent])

        system_params['rigidWater'] = False  # not currently supported (because I'm lazy)
        system_params['constraints'] = None
        if self.mol.integrator is not None:
            if self.mol.integrator.params.get('constrain_water', False):
                system_params['rigidWater'] = True
            if self.mol.integrator.params.get('constrain_hbonds', False):
                system_params['constraints'] = app.HBonds

        # Deal with h-bonds listed in molecular constraints
        for constraint in self.mol.constraints:
            if constraint.desc == 'hbonds':
                system_params['constraints'] = app.HBonds
                break

        return system_params

    def _get_platform(self):
        from simtk import openmm

        preference = ['CUDA', 'OpenCL', 'CPU', 'Reference']
        from_lower = {x.lower(): x for x in preference}

        properties = {}

        if self.params.compute_platform.lower() == 'auto':
            for platname in preference:
                try:
                    platform = openmm.Platform.getPlatformByName(platname)
                except Exception:  # it just throws "Exception" unfortunately
                    continue
                else:
                    self.params.compute_platform = platname.lower()
                    break

            raise moldesign.NotSupportedError("Likely OpenMM installation error. "
                                              "none of the expected platforms were found: "
                                              + ', '.join(preference))
        else:
            platform = openmm.Platform.getPlatformByName(from_lower[self.params.compute_platform])

        if self.params.compute_platform == 'cpu' and self.params.num_cpus > 0:
            # need to use native_strs here or the swig interface gets confused
            properties[native_str('Threads')] = native_str(self.params.num_cpus)

        return platform, properties
