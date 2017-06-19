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

import itertools
import imp

import numpy as np

import pyccc
import moldesign as mdt
from ..utils import from_filepath
from .. import units as u
from .. import compute
from ..utils import exports


try:
    imp.find_module('simtk')
except (ImportError, OSError) as exc:
    print('OpenMM could not be imported; using remote docker container')
    force_remote = True
else:
    force_remote = False


class OpenMMPickleMixin(object):
    def __getstate__(self):
        mystate = self.__dict__.copy()
        if 'sim' in mystate:
            assert 'sim_args' not in mystate
            sim = mystate.pop('sim')
            mystate['sim_args'] = (sim.topology, sim.system, sim.integrator)
        return mystate

    def __setstate__(self, state):
        from simtk.openmm import app
        if 'sim_args' in state:
            assert 'sim' not in state
            args = state.pop('sim_args')
            state['sim'] = app.Simulation(*args)
        self.__dict__.update(state)


# This is a factory for the MdtReporter class. It's here so that we don't have to import
# simtk.openmm.app at the module level
def MdtReporter(mol, report_interval):
    from simtk.openmm.app import StateDataReporter

    class MdtReporter(StateDataReporter):
        """
        We'll use this class to capture all the information we need about a trajectory
        It's pretty basic - the assumption is that there will be more processing on the client side
        """

        LEN = 30
        def __init__(self, mol, report_interval):
            self.mol = mol
            self.report_interval = report_interval
            self.trajectory = mdt.Trajectory(mol)
            self.annotation = None
            self._row_format = ("{:<%d}" % 10) + 3*("{:>%d}" % self.LEN)
            self._printed_header = False
            self.last_report_time = None

        def report_from_mol(self, **kwargs):
            self.mol.calculate()
            if self.annotation is not None:
                kwargs.setdefault('annotation', self.annotation)
            self.report(self.mol.energy_model.sim,
                        self.mol.energy_model.sim.context.getState(getEnergy=True,
                                                                   getForces=True,
                                                                   getPositions=True,
                                                                   getVelocities=True))

        def report(self, simulation, state):
            """ Callback for dynamics after the specified interval

            Args:
               simulation (simtk.openmm.app.Simulation): simulation to report on
               state (simtk.openmm.State): state of the simulation
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
                print(self._row_format.format(timeheader, peheader, keheader, temperatureheader))
                self._printed_header = True
            ke = mdt.helpers.kinetic_energy(report['momenta'], self.mol.dim_masses)
            t = (2.0 * ke) / (u.k_b * self.mol.dynamic_dof)
            print(self._row_format.format(report['time'].defunits_value(),
                                          report['potential_energy'].defunits_value(),
                                          ke.defunits_value(),
                                          t.defunits_value()))
            self.last_report_time = self.mol.time

            self.trajectory.new_frame(properties=report)

        def describeNextReport(self, simulation):
            """
            Returns:
                tuple: A five element tuple.  The first element is the number of steps
                until the next report.  The remaining elements specify whether
                that report will require positions, velocities, forces, and
                energies respectively.
            """
            steps = self.report_interval - simulation.currentStep % self.report_interval
            return (steps, True, True, True, True)

    return MdtReporter(mol, report_interval)


PINT_NAMES = {'mole': u.avogadro,
              'degree': u.degrees,
              'radian': u.radians,
              'elementary charge': u.q_e}


@exports
def simtk2pint(quantity, flat=False):
    """ Converts a quantity from the simtk unit system to the internal unit system

    Args:
        quantity (simtk.unit.quantity.Quantity): quantity to convert
        flat (bool): if True, flatten 3xN arrays to 3N

    Returns:
        mdt.units.MdtQuantity: converted to MDT unit system
    """
    from simtk import unit as stku

    mag = np.array(quantity._value)

    if quantity.unit == stku.radian:
        return mag * u.radians
    if quantity.unit == stku.degree:
        return mag * u.degrees

    for dim, exp in itertools.chain(quantity.unit.iter_scaled_units(),
                                    quantity.unit.iter_top_base_units()):
        if dim.name in PINT_NAMES:
            pintunit = PINT_NAMES[dim.name]
        else:
            pintunit = u.ureg.parse_expression(dim.name)
        mag = mag * (pintunit**exp)
        if flat:
            mag = np.reshape(mag, (np.product(mag.shape),))
    return u.default.convert(mag)


@exports
def pint2simtk(quantity):
    """ Converts a quantity from the pint to simtk unit system.
    Note SimTK has a less extensive collection that pint. May need to have pint convert
    to SI first
    """
    from simtk import unit as stku

    SIMTK_NAMES = {'ang': stku.angstrom,
                   'fs': stku.femtosecond,
                   'nm': stku.nanometer,
                   'ps': stku.picosecond}

    newvar = quantity._magnitude
    for dim, exp in quantity._units.items():
        if dim in SIMTK_NAMES:
            stkunit = SIMTK_NAMES[dim]
        else:
            stkunit = getattr(stku, dim)
        newvar = newvar * stkunit ** exp
    return newvar


@compute.runsremotely(enable=force_remote)
def _amber_to_mol(prmtop_file, inpcrd_file):
    """ Convert an amber prmtop and inpcrd file to an MDT molecule

    Args:
        prmtop_file (file-like): topology file in amber prmtop format
        inpcrd_file (file-like): coordinate file in amber crd format

    Returns:
        moldesign.Molecule: Molecule parsed from amber output
    """
    from simtk.openmm import app

    prmtop = from_filepath(app.AmberPrmtopFile, prmtop_file)
    inpcrd = from_filepath(app.AmberInpcrdFile, inpcrd_file)

    mol = topology_to_mol(prmtop.topology,
                          positions=inpcrd.positions,
                          velocities=inpcrd.velocities)
    return mol


if force_remote:
    def amber_to_mol(prmtop_file, inpcrd_file):
        if not isinstance(prmtop_file, pyccc.FileContainer):
            prmtop_file = pyccc.LocalFile(prmtop_file)
        if not isinstance(inpcrd_file, pyccc.FileContainer):
            inpcrd_file = pyccc.LocalFile(inpcrd_file)
        return _amber_to_mol(prmtop_file, inpcrd_file)
else:
    amber_to_mol = _amber_to_mol


exports(amber_to_mol)


@exports
def topology_to_mol(topo, name=None, positions=None, velocities=None, assign_bond_orders=True):
    """ Convert an OpenMM topology object into an MDT molecule.

    Args:
        topo (simtk.openmm.app.topology.Topology): topology to convert
        name (str): name to assign to molecule
        positions (list): simtk list of atomic positions
        velocities (list): simtk list of atomic velocities
        assign_bond_orders (bool): assign bond orders from templates (simtk topologies
             do not store bond orders)

    """
    from simtk import unit as stku

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
        na1, na2 = atommap[a1], atommap[a2]
        if na1 not in bonds:
            bonds[na1] = {}
        if na2 not in bonds:
            bonds[na2] = {}
        bonds[na1][na2] = 1
        bonds[na2][na1] = 1

    if name is None:
        name = 'Unnamed molecule from OpenMM'

    newmol = mdt.Molecule(newatoms, bond_graph=bonds, name=name)

    if assign_bond_orders:
        for residue in newmol.residues:
            try:
                residue.assign_template_bonds()
            except (KeyError, ValueError):
                pass

    return newmol


@exports
def mol_to_topology(mol):
    """ Create an openmm topology object from an MDT molecule

    Args:
        mol (moldesign.Molecule): molecule to copy topology from

    Returns:
        simtk.openmm.app.Topology: topology of the molecule
    """
    from simtk.openmm import app

    top = app.Topology()
    chainmap = {chain: top.addChain(chain.name) for chain in mol.chains}
    resmap = {res: top.addResidue(res.resname, chainmap[res.chain], str(res.pdbindex))
              for res in mol.residues}
    atommap = {atom: top.addAtom(atom.name,
                                 app.Element.getBySymbol(atom.element),
                                 resmap[atom.residue],
                                 id=str(atom.pdbindex))
               for atom in mol.atoms}
    for bond in mol.bonds:
        top.addBond(atommap[bond.a1], atommap[bond.a2])

    return top


@exports
def mol_to_modeller(mol):
    from simtk.openmm import app
    if mol.is_small_molecule:
        if not mol.residues[0].resname:
            mol.residues[0].resname = 'UNL'
            mol.residues[0].pdbindex = 1
        if not mol.chains[0].pdbname:
            mol.chains[0].pdbname = 'A'

    return app.Modeller(mol_to_topology(mol), pint2simtk(mol.positions))

