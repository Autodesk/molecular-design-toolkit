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

import numpy as np

import pyccc
from moldesign.utils import from_filepath

try:
    import simtk.openmm as mm
    from simtk import unit as stku
    from simtk.openmm import app
except ImportError:
    mm = stku = app = None
    force_remote = True
else:
    force_remote = False

from pyccc import LocalFile

import moldesign as mdt
from moldesign import units as u

from moldesign import compute

from moldesign.molecules import Trajectory, Molecule


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


# This class needs special treatment because it inherits from an
# openmm class, but OpenMM may not be present in the user's environment
if force_remote:
    MdtReporter = None
else:
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

            # TODO: make sure openmm masses are the same as MDT masses
            report = dict(
                positions=simtk2pint(state.getPositions()),
                momenta=simtk2pint(state.getVelocities())*self.mol.dim_masses,
                forces=simtk2pint(state.getForces()),
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
            ke = mdt.helpers.kinetic_energy(report['momenta'], self.mol.dim_masses)
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


@compute.runsremotely(enable=force_remote)
def _amber_to_mol(prmtop_file, inpcrd_file):
    prmtop = from_filepath(app.AmberPrmtopFile, prmtop_file)
    inpcrd = from_filepath(app.AmberInpcrdFile, inpcrd_file)

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
