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
import itertools
import tempfile
from cStringIO import StringIO

import numpy as np

import pyccc

try:
    import simtk.openmm as mm
    from simtk import unit as stku
    from simtk.openmm import app
except ImportError:
    force_remote = True
else:
    force_remote = False

from pyccc import LocalFile

import moldesign as mdt
from moldesign import units as u

from moldesign import forcefields as ff
from moldesign import compute

from moldesign.integrators.base import IntegratorBase, LangevinBase
from moldesign.models.base import MMBase
from moldesign.structure import Trajectory, MolecularProperties, Molecule


class OpenMMPickleMixin(object):
    def __getstate__(self):
        mystate = self.__dict__.copy()
        if 'sim' in mystate:
            assert 'sim_args' not in mystate
            sim = mystate.pop('sim')
            mystate['sim_args'] = (sim.topology, sim.system, sim.integrator)
        return mystate

    def __setstate__(self, state):
        if 'sim_args' in state:
            assert 'sim' not in state
            args = state.pop('sim_args')
            state['sim'] = app.Simulation(*args)
        self.__dict__.update(state)


class OpenMMBaseIntegrator(IntegratorBase, OpenMMPickleMixin):
    _openmm_compatible = True

    def prep(self):
        assert self.mol.energy_model._openmm_compatible
        if self._prepped and self.model is self.mol.energy_model and self.model._prepped: return
        self.model = self.mol.energy_model
        self.model.prep()
        self.sim = self.model.sim
        self.reporter = self._attach_reporters()
        self._prepped = True

    def run(self, run_for, wait=False):
        # TODO: like model.minimize, this is a hacky wrapper that we need to replace with
        # something more generalizable
        try:
            traj = self._run(run_for)
        except pyccc.ProgramFailure:
            raise pyccc.ProgramFailure('OpenMM crashed silently. Please examine the output. '
                                       'This may be due to large forces from, for example, '
                                       'an insufficiently minimized starting geometry.')
        if force_remote or (not wait):
            self.mol.energy_model._sync_remote(traj.mol)
            traj.mol = self.mol
        return traj

    @compute.runsremotely(remote=force_remote, is_imethod=True)
    def _run(self, run_for):
        self.prep()
        nsteps = self.time_to_steps(run_for, self.params.timestep)
        self.mol.time = 0.0 * u.default.time
        self.model._set_openmm_state()
        self.reporter.annotation = 'Langevin dynamics @ %s' % self.params.temperature
        self.reporter.report_from_mol()
        self.sim.step(nsteps)
        if self.reporter.last_report_time != self.mol.time:
            self.reporter.report_from_mol()
        self.reporter.report_from_mol()
        self.model._sync_to_openmm()
        return self.reporter.trajectory

    def _attach_reporters(self):
        """
        Make sure the simulation has reporters for this run
        :return:
        """
        report_interval = self.time_to_steps(self.params.frame_interval,
                                             self.params.timestep)
        reporter = MdtReporter(self.mol, report_interval)
        self.sim.reporters = [reporter]
        return reporter


class OpenMMVerlet(OpenMMBaseIntegrator):
    def get_openmm_integrator(self):
        integrator = mm.VerletIntegrator(pint2simtk(self.params.timestep))
        return integrator


class OpenMMLangevin(LangevinBase, OpenMMBaseIntegrator):
    def get_openmm_integrator(self):
        integrator = mm.LangevinIntegrator(
            pint2simtk(self.params.temperature),
            pint2simtk(self.params.collision_rate),
            pint2simtk(self.params.timestep))
        return integrator


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
                                     'bond (%s)' % mdt.Bond(ai, aj))
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


### This class needs special treatment because it inherits from an
### openmm class, but OpenMM may not be present in the user's environment
if not force_remote:
    class MdtReporter(app.StateDataReporter):
        """
        We'll use this class to capture all the information we need about a trajectory
        It's pretty basic - the assumption is that there will be more processing on the client side
        """

        LEN = 30
        def __init__(self, mol, report_interval):
            self.mol = mol
            self.report_interval = report_interval
            self.trajectory = Trajectory(mol)
            self.annotation = None
            self._row_format = ("{:<%d}" % 10) + 3*("{:>%d}" % self.LEN)
            self._printed_header = False
            self.last_report_time = None

        def report_from_mol(self, **kwargs):
            self.mol.calculate()
            if self.annotation is not None:
                kwargs.setdefault('annotation',self.annotation)
            self.trajectory.new_frame(**kwargs)

        def report(self, simulation, state):
            """
            Callback for dynamics after the specified interval
            :type simulation: simtk.openmm.app.Simulation
            :type state: simtk.openmm.State
            """
            report = dict(
                positions=simtk2pint(state.getPositions(), flat=True),
                momenta=simtk2pint(state.getVelocities(), flat=True)*self.mol.dim_masses,
                forces=simtk2pint(state.getForces(), flat=True),
                time=simtk2pint(state.getTime()),
                vectors=simtk2pint(state.getPeriodicBoxVectors()),
                potential_energy=simtk2pint(state.getPotentialEnergy()))
            if self.annotation is not None: report['annotation'] = self.annotation

            if not self._printed_header:
                timeheader = 'time / {units}'.format(units=u.default.time)
                peheader = 'potential / {units}'.format(units=u.default.energy)
                keheader = 'kinetic / {units}'.format(units=u.default.energy)
                temperatureheader = 'T / {units}'.format(units=u.default.temperature)
                print self._row_format.format(timeheader, peheader, keheader, temperatureheader)
                self._printed_header = True
            ke = (report['momenta'] / (2*self.mol.dim_masses)).dot(report['momenta'])
            t = (2.0 * ke) / (u.k_b * self.mol.dynamic_dof)
            print self._row_format.format(report['time'].defunits_value(),
                                          report['potential_energy'].defunits_value(),
                                          ke.defunits_value(),
                                          t.defunits_value())
            self.last_report_time = self.mol.time


            self.trajectory.new_frame(properties=report)

        def describeNextReport(self, simulation):
            """
            :type simulation: simtk.openmm.app.Simulation
            :return: report_description
                A five element tuple.  The first element is the number of steps
                until the next report.  The remaining elements specify whether
                that report will require positions, velocities, forces, and
                energies respectively.
            :return type: tuple
            """
            steps = self.report_interval - simulation.currentStep % self.report_interval
            return (steps, True, True, True, True)


PINT_NAMES = {'mole': u.avogadro,
              'degree': u.degrees,
              'radian': u.radians,
              'elementary charge': u.q_e}


def simtk2pint(quantity, flat=False):
    """
    Converts a quantity from the simtk unit system to a quantity from the pint unit system
    :param quantity:
    :param flat: if True, flatten 3xN arrays to 3N
    """
    mag = quantity._value

    if quantity.unit == stku.radian:
        return mag * u.radians
    if quantity.unit == stku.degree:
        return mag * u.degrees

    if hasattr(mag, '__getslice__'): mag = np.array(mag[:])
    for dim, exp in itertools.chain(quantity.unit.iter_scaled_units(),
                                    quantity.unit.iter_top_base_units()):
        if dim.name in PINT_NAMES:
            pintunit = PINT_NAMES[dim.name]
        else:
            pintunit = u.ureg.parse_expression(dim.name)
        mag = mag * (pintunit**exp)
        if flat: mag = np.reshape(mag, (np.product(mag.shape),))
    return u.default.convert(mag)


def pint2simtk(quantity):
    """ Converts a quantity from the pint to simtk unit system.
    Note SimTK appears limited, esp for energy units. May need to have pint convert
    to SI first
    """
    SIMTK_NAMES = {'ang': stku.angstrom,
                   'fs': stku.femtosecond,
                   'nm': stku.nanometer,
                   'ps': stku.picosecond}

    newvar = quantity._magnitude
    for dim, exp in quantity._units.iteritems():
        if dim in SIMTK_NAMES:
            stkunit = SIMTK_NAMES[dim]
        else:
            stkunit = getattr(stku, dim)
        newvar = newvar * stkunit ** exp
    return newvar


def from_filepath(filelike, func):
    """Run on func on a temporary *path* assigned to filelike"""
    if type(filelike) == str:
        return func(filelike)
    else:
        with tempfile.NamedTemporaryFile() as outfile:
            outfile.write(filelike.read())
            outfile.flush()
            result = func(outfile.name)
        return result


@compute.runsremotely(remote=force_remote)
def _amber_to_mol(prmtop_file, inpcrd_file):
    prmtop = from_filepath(prmtop_file, app.AmberPrmtopFile)
    inpcrd = from_filepath(inpcrd_file, app.AmberInpcrdFile)

    mol = topology_to_mol(prmtop.topology,
                          positions=inpcrd.positions,
                          velocities=inpcrd.velocities)
    return mol


if force_remote:
    def amber_to_mol(prmtop_file, inpcrd_file):
        if not isinstance(prmtop_file, pyccc.FileContainer):
            prmtop_file = LocalFile(prmtop_file)
        if not isinstance(inpcrd_file, pyccc.FileContainer):
            inpcrd_file = LocalFile(inpcrd_file)
        return _amber_to_mol(prmtop_file, inpcrd_file)
else:
    amber_to_mol = _amber_to_mol


def list_openmmplatforms():
    return [mm.Platform.getPlatform(ip).getName()
            for ip in xrange(mm.Platform.getNumPlatforms())]


def topology_to_mol(topo, name=None, positions=None, velocities=None):
    """
    :type topo: simtk.openmm.app.topology.Topology
    :return:
    """
    # Atoms
    atommap = {}
    newatoms = []
    masses = u.amu*[atom.element.mass.value_in_unit(stku.amu) for atom in topo.atoms()]
    for atom,mass in zip(topo.atoms(), masses):
        newatom = mdt.Atom(atnum=atom.element.atomic_number,
                           name=atom.name,
                           mass=mass)
        atommap[atom] = newatom
        newatoms.append(newatom)

    # Coordinates
    if positions is not None:
        poslist = np.array([p.value_in_unit(stku.nanometer) for p in positions]) * u.nm
        poslist.ito(u.default.length)
        for newatom, position in zip(newatoms, poslist):
            newatom.position = position
    if velocities is not None:
        velolist = np.array([v.value_in_unit(stku.nanometer/stku.femtosecond) for v in velocities]) * u.nm/u.fs
        velolist = u.default.convert(velolist)
        for newatom, velocity in zip(newatoms, velolist):
            newatom.momentum = newatom.mass * simtk2pint(velocity)

    # Biounits
    chains = {}
    for chain in topo.chains():
        if chain.id not in chains:
            chains[chain.id] = mdt.Chain(name=chain.id, index=chain.index)
        newchain = chains[chain.id]
        for residue in chain.residues():
            newresidue = mdt.Residue(name='%s%d' % (residue.name,
                                                         residue.index),
                                          chain=newchain,
                                          pdbindex=int(residue.id),
                                          pdbname=residue.name)
            newchain.add(newresidue)
            for atom in residue.atoms():
                newatom = atommap[atom]
                newatom.residue = newresidue
                newatom.chain = newchain
                newresidue.add(newatom)

    # Bonds
    bonds = {}
    for bond in topo.bonds():
        a1, a2 = bond
        na1, na2 = atommap[a1],atommap[a2]
        if na1 not in bonds:
            bonds[na1] = {}
        if na2 not in bonds:
            bonds[na2] = {}
        bonds[na1][na2] = -1
        bonds[na2][na1] = -1

    if name is None: name = 'Unnamed molecule from OpenMM'
    newmol = Molecule(newatoms, bond_graph=bonds, name=name)
    return newmol


def mol_to_toplogy(mol):
    raise NotImplementedError


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

