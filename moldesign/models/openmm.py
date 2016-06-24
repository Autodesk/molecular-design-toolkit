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

from cStringIO import StringIO

import moldesign.structure
from moldesign import compute
from moldesign import forcefields as ff
from moldesign.structure import Trajectory, MolecularProperties

from moldesign.interfaces.openmm import OpenMMPickleMixin, force_remote, simtk2pint, \
    pint2simtk, from_filepath, app, stku, mm

from .base import MMBase


def exports(o):
    __all__.append(o.__name__)
    return o
__all__ = []


@exports
class OpenMMPotential(MMBase, OpenMMPickleMixin):
    """Creates an OpenMM "context" to drive energy calculations.
    Note that, while a dummy integrator is assigned, a different context will
    be created for any MD calculations.

    :ivar sim: openmm simulation object
    :type sim: simtk.openmm.app.Simulation
    """
    # NEWFEATURE: need to set/get platform (and properties, e.g. number of threads)
    DEFAULT_PROPERTIES = ['potential_energy', 'forces']
    _openmm_compatible = True

    def __init__(self, **kwargs):
        super(OpenMMPotential, self).__init__(**kwargs)
        self.sim = None

    def get_openmm_simulation(self):
        if force_remote:
            raise ImportError("Can't create an OpenMM object on this machine - OpenMM not "
                              "installed")
        else:
            if not self._prepped: self.prep()
            return self.sim

    @compute.runsremotely(remote=force_remote, is_imethod=True)
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
                                    potential_energy=simtk2pint(state.getPotentialEnergy()),
                                    forces=simtk2pint(state.getForces(), flat=True))
        props['forces'].shape = (self.mol.ndims,)
        return props

    def prep(self, force=False):
        """
        Drive the construction of the openmm simulation
        This will rebuild this OpenMM simulation if: A) it's not built yet, or B) there's a new integrator
        """

        # TODO: automatically set _prepped to false if the model or integration parameters change
        if not force:
            if self._prepped and self._prep_integrator == self.mol.integrator: return
        self._create_system()

        if self.mol.integrator is None:
            integrator = self._make_dummy_integrator()
        else:
            integrator = self.mol.integrator.get_openmm_integrator()

        self._set_constraints()
        self.sim = app.Simulation(self.mm_topology, self.mm_system, integrator)
        self._prepped = True
        print 'Created OpenMM kernel (Platform: %s)' % self.sim.context.getPlatform().getName()
        self._prep_integrator = self.mol.integrator

    def minimize(self, **kwargs):
        traj = self._minimize(**kwargs)

        if force_remote or (not kwargs.get('wait', False)):
            self._sync_remote(traj.mol)
            traj.mol = self.mol

        return traj

    def _sync_remote(self, mol):
        # TODO: this is a hack to update the object after a minimization
        #        We need a better pattern for this, ideally one that doesn't
        #        require an explicit wrapper like this - we shouldn't have to copy
        #        the properties over manually
        self.mol.positions = mol.positions
        self.mol.momenta = mol.momenta
        self.mol.properties = mol.properties
        self.mol.time = mol.time


    @compute.runsremotely(remote=force_remote, is_imethod=True)
    def _minimize(self, nsteps=500,
                 force_tolerance=None,
                 frame_interval=None,
                 wait=True):
        """ Run an OpenMM minimization.
        Note:
            We're not able to support `frame_interval` for this method;
                see https://github.com/pandegroup/openmm/issues/1155
        Args:
           nsteps(int): maximum number of steps
           force_tolerance (moldesign.units.MdtQuantity): RMSD tolerance for convergence [energy/length]
        """

        # NEWFEATURE: write/find an openmm "integrator" to do this minimization
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

    @compute.runsremotely(remote=force_remote, is_imethod=True)
    def get_forcefield(self):
        """ Get the force field parameters for this molecule.

        Note:
            The returned object is for introspection only; it can't be used (yet) to modify the energy model, and
                it doesn't include the entire force field (most notably missing are periodic replicate forces,
                1-4 nonbonded attenuation, and implicit solvent effects)

        Returns:
            moldesign.forcefield.ForceField
        """
        if self.mdtforcefield is None:
            if self._prepped:
                self.mdtforcefield = _get_system_params(self.sim.system, self.mol)

            elif not hasattr(self, 'mm_system'):
                self._create_system()

            self.mdtforcefield = _get_system_params(self.mm_system, self.mol)

        return self.mdtforcefield


    #################################################
    # "Private" methods for managing OpenMM are below
    def _set_constraints(self):
        system = self.mm_system
        fixed_atoms = set()

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
                system.addConstraint(constraint.a1,
                                     constraint.a2,
                                     pint2simtk(constraint.value))

            else:
                raise ValueError('OpenMM interface does not support "%s" constraints.' % constraint.desc)

        # Workaround for OpenMM issue: can't have an atom that's both fixed *and* has a distance constraint.
        # If both atoms in the distance constraint are also fixed, then we can just remove the constraint
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
                elif (ai in fixed_atoms) or (aj in fixed_atoms):  #only one is fixed
                    raise ValueError('In OpenMM, fixed atoms cannot be part of a constrained '
                                     'bond (%s)'%moldesign.structure.bonds.Bond(ai, aj))
                else:
                    ic += 1

    @staticmethod
    def _make_dummy_integrator():
        return mm.VerletIntegrator(2.0 * stku.femtoseconds)

    def _create_system(self):
        # Parse the stored PRMTOP file if it's available
        if ('amber_params' in self.mol.ff) and not hasattr(self, 'mm_prmtop'):
            print 'Parsing stored PRMTOP file: %s' % self.mol.ff.amber_params.prmtop
            self.mm_prmtop = from_filepath(self.mol.ff.amber_params.prmtop, app.AmberPrmtopFile)

        # Create the OpenMM system
        system_params = self._get_system_params()
        if hasattr(self, 'mm_prmtop'):
            self.mm_system = self.mm_prmtop.createSystem(**system_params)
            self.mm_topology = self.mm_prmtop.topology
            print 'Created force field using embedded prmtop file'
        else:
            self.mm_pdb = from_filepath(StringIO(self.mol.write(format='pdb')), app.PDBFile)
            self.mm_topology = self.mm_pdb.topology
            ffs = self._get_forcefield_args()
            self.mm_forcefield = app.ForceField(*self._get_forcefield_args())
            print 'Created force field using OpenMM templates: %s' % str(ffs)
            self.mm_system.createSystem(self.mm_topology, self.mol.name,
                                        **system_params)

    def _get_forcefield_args(self):
        ff = self.params.forcefield[:2]
        if ff[:5] == 'amber':
            return ('%s.xml' % ff, 'tip3p.xml')
        else:
            raise NotImplementedError()

    def _set_openmm_state(self):
        self.sim.context.setPositions(pint2simtk(self.mol.atoms.position))
        self.sim.context.setVelocities(pint2simtk(self.mol.atoms.velocity))

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
                    atom.position = simtk2pint(positions[iatom])
                if velocities != False:
                    atom.momentum = simtk2pint(velocities[iatom]) * atom.mass

        if time is True:
            if state is None:
                state = self.sim.context.getState()
            time = simtk2pint(state.getTime())
        elif time:
            time = simtk2pint(time)
        if time:
            self.time = time

    def _get_system_params(self):
        """
        Translates the spec from MMBase into system parameter keywords
        for createSystem
        """
        # need cmm motion
        nonbonded_names = {'nocutoff': app.NoCutoff,
                           'ewald': app.Ewald,
                           'pme': app.PME,
                           'cutoff': app.CutoffPeriodic if self.params.periodic else app.CutoffNonPeriodic}
        implicit_solvent_names = {'obc': app.OBC2}

        system_params = dict(nonbondedMethod=nonbonded_names[self.params.nonbonded],
                             nonbondedCutoff=pint2simtk(self.params.cutoff),
                             implicitSolvent=implicit_solvent_names[self.params.implicit_solvent])

        if self.params.get('constrain_water', False):
            system_params['rigidWater'] = True
        if self.params.get('constrain_hbonds', False):
            system_params['constraints'] = app.HBonds

        return system_params


def list_openmmplatforms():
    return [mm.Platform.getPlatform(ip).getName()
            for ip in xrange(mm.Platform.getNumPlatforms())]



def _get_system_params(system, mol):
    # TODO: 1-4 bond rules
    # TODO: constraints
    forces = system.getForces()
    bonds, angles, dihedrals = [], [], []
    charges = {}
    ljparameters = {}

    for f in forces:
        if type(f) == mm.HarmonicBondForce:
            for ibond in xrange(f.getNumBonds()):
                i1, i2, d0, k = f.getBondParameters(ibond)
                bond = ff.HarmonicBondTerm(mol.atoms[i1], mol.atoms[i2],
                                           simtk2pint(k),
                                           simtk2pint(d0))
                bonds.append(bond)

        elif type(f) == mm.HarmonicAngleForce:
            for iangle in xrange(f.getNumAngles()):
                i1, i2, i3, t0, k = f.getAngleParameters(iangle)
                angle = ff.HarmonicAngleTerm(mol.atoms[i1], mol.atoms[i2], mol.atoms[i3],
                                             simtk2pint(k),
                                             simtk2pint(t0))
                angles.append(angle)

        elif type(f) == mm.PeriodicTorsionForce:
            for itorsion in xrange(f.getNumTorsions()):
                i1, i2, i3, i4, n, t0, v_n = f.getTorsionParameters(itorsion)
                torsion = ff.PeriodicTorsionTerm(mol.atoms[i1], mol.atoms[i2], mol.atoms[i3], mol.atoms[i4],
                                                 n, simtk2pint(v_n),
                                                 simtk2pint(t0))
                dihedrals.append(torsion)

        elif type(f) == mm.NonbondedForce:
            for i, atom in enumerate(mol.atoms):
                q, sigma, epsilon = f.getParticleParameters(i)
                charges[atom] = simtk2pint(q)
                ljparameters[atom] = ff.LennardJonesSigmaEps(atom,
                                                             simtk2pint(sigma),
                                                             simtk2pint(epsilon))

        elif type(f) == mm.CMMotionRemover:
            continue

        elif type(f) == mm.GBSAOBCForce:
            continue

    ffparams = ff.FFParameters(bonds, angles, dihedrals, charges, ljparameters)
    return ffparams

