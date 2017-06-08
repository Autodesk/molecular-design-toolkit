"""
This module contains various helper functions used by MDT internally.
"""
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

import collections

import numpy as np

from moldesign import units as u


class VolumetricGrid(object):
    """
    Helper object for preparing gaussian CUBE files
    """
    UNITS = u.angstrom
    def __init__(self, positions, padding=2.5*u.angstrom, npoints=25):
        mins = positions.min(axis=0) - padding
        maxes = positions.max(axis=0) + padding
        self.npoints = npoints
        self.xr = (mins[0], maxes[0])
        self.yr = (mins[1], maxes[1])
        self.zr = (mins[2], maxes[2])
        self._origin = mins.value_in(self.UNITS)
        self.dx = (self.xr[1] - self.xr[0]).value_in(self.UNITS) / (float(npoints) - 1)
        self.dy = (self.yr[1] - self.yr[0]).value_in(self.UNITS) / (float(npoints) - 1)
        self.dz = (self.zr[1] - self.zr[0]).value_in(self.UNITS) / (float(npoints) - 1)
        self.fxyz = None

    def xyzlist(self):
        stride = self.npoints * 1j
        grids = np.mgrid[self.xr[0]:self.xr[1]:stride,
                self.yr[0]:self.yr[1]:stride,
                self.zr[0]:self.zr[1]:stride]
        return grids * self.UNITS

    def origin(self):
        return tuple(self._origin)


def get_all_atoms(*objects):
    """ Given Atoms, AtomContainers, lists of Atoms, and lists of AtomContainers,
    return a flat list of all atoms contained therein.

    A given atom is only returned once, even if it's found more than once.

    Args:
        *objects (moldesign.Atom OR moldesign.AtomContainer OR List[moldesign.Atom] OR
            List[moldesign.AtomContainer]): objects to take atoms from

    """
    from moldesign import molecules

    atoms = collections.OrderedDict()

    for obj in objects:
        if isinstance(obj, molecules.Atom):
            atoms[obj] = None
        elif hasattr(obj, 'atoms'):
            atoms.update((x,None) for x in obj.atoms)
        else:
            for item in obj:
                if isinstance(item, molecules.Atom):
                    atoms[item] = None
                elif hasattr(item, 'atoms'):
                    atoms.update((x, None) for x in item.atoms)

    return molecules.AtomList(iter(atoms.keys()))


def kinetic_energy(momenta, masses):
    return 0.5 * (momenta*momenta/masses).sum()


def kinetic_temperature(ke, dof):
    from moldesign.units import k_b
    t = (2.0*ke)/(k_b*dof)
    return t.defunits()


def atom_name_check(mol, force=False):
    """ Makes sure atom names are unique in each residue.

    If atoms names aren't unqiue:
      - if the names are just the names of the elements, rename them
      - else print a warning
    """
    badres = []
    for residue in mol.residues:
        names = set(atom.name for atom in residue.atoms)
        if len(names) != residue.num_atoms:
            # atom names aren't unique, check if we can change them
            for atom in residue.atoms:
                if atom.name.lower() != atom.symbol.lower():
                    badres.append(residue)
                    if not force:
                        break
            else:  # rename the atoms
                atomnums = {}
                for atom in residue.atoms:
                    atom.name = atom.symbol + str(atomnums.setdefault(atom.symbol, 0))
                    atomnums[atom.symbol] += 1

    if badres:
        print('WARNING: residues do not have uniquely named atoms: %s' % badres)


def restore_topology(mol, topo):
    """ Restores chain IDs and residue indices (these are stripped by some methods)

    Args:
        mol (mdt.Molecule): molecule to restore topology to
        topo (mdt.Molecule): reference topology

    Returns:
        mdt.Molecule: a copy of ``mol`` with a restored topology
    """
    import moldesign as mdt

    assert mol.num_residues == topo.num_residues
    assert mol.num_chains == 1

    chain_map = {}
    for chain in topo.chains:
        chain_map[chain] = mdt.Chain(name=chain.name)

    for res, refres in zip(mol.residues, topo.residues):
        if refres.resname in ('HID', 'HIE', 'HIP'):
            rname = 'HIS'
        else:
            rname = refres.resname
        assert res.resname == rname
        res.pdbindex = refres.pdbindex
        res.name = refres.name
        res.chain = chain_map[refres.chain]
        for atom in res.atoms:
            atom.chain = res.chain

    return mdt.Molecule(mol.atoms)
