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

import imp
import io

import numpy as np

import moldesign as mdt
from .. import units as u
from .. import compute
from ..utils import exports
from . import openmm as opm

try:
    imp.find_module('pdbfixer')
except (ImportError, OSError) as exc:
    print('PDBFixer could not be imported; using remote docker container')
    force_remote = True
else:
    force_remote = False


@exports
def mol_to_fixer(mol):
    import pdbfixer
    fixer = pdbfixer.PDBFixer(
            pdbfile=io.BytesIO(mol.write(format='pdb')))
    return fixer


@exports
def fixer_to_mol(f):
    return opm.topology_to_mol(f.topology, positions=f.positions)


@compute.runsremotely(enable=force_remote)
def add_hydrogen(mol, pH=7.4):
    fixer = mol_to_fixer(mol)
    fixer.addMissingHydrogens(pH)
    return fixer_to_mol(fixer)


@compute.runsremotely(enable=force_remote)
def mutate(mol, mutationmap):
    fixer = mol_to_fixer(mol)
    chain_mutations = {}
    for res in mutationmap:
        chain_mutations.setdefault(res.chain.pdbname, {})[res] = mutationmap[res]

    for chainid, mutations in chain_mutations.items():
        mutstrings = ['%s-%d-%s' % (res.resname, res.pdbindex, newname)
                      for res, newname in mutations.items()]
        print('Applying mutations to chain %s: %s' % (chainid, ', '.join(mutstrings)))
        fixer.applyMutations(mutations, chainid)
    return fixer_to_mol(fixer)


@compute.runsremotely(enable=force_remote)
def get_missing_residues(mol):
    fixer = mol_to_fixer(mol)
    fixer.findMissingResidues()
    fixerchains = list(fixer.topology.chains)

    missing = list()
    for (chainidx, insertionpoint), reslist in fixer.missingResidues.items():
        chainid = fixerchains[chainidx].id
        for ires, resname in enumerate(reslist):
            missing.append(mdt.helpers.MissingResidue(chainid, resname, insertionpoint + ires))

    return missing


@compute.runsremotely(enable=force_remote)
def add_water(mol, min_box_size=None, padding=None,
              ion_concentration=None, neutralize=True,
              positive_ion='Na+', negative_ion='Cl-'):
    """ Solvate a molecule in a water box with optional ions

    Args:
        mol (moldesign.Molecule): solute molecule
        min_box_size (u.Scalar[length] or u.Vector[length]): size of the water box - either
           a vector of x,y,z dimensions, or just a uniform cube length. Either this or
           ``padding`` (or both) must be passed
        padding (u.Scalar[length]): distance to edge of water box from the solute in each dimension
        neutralize (bool): add ions to neutralize solute charge (in
            addition to specified ion concentration)
        positive_ion (str): type of positive ions to add, if needed. Allowed values
            (from OpenMM modeller) are Cs, K, Li, Na (the default) and Rb
        negative_ion (str): type of negative ions to add, if needed. Allowed values
            (from OpenMM modeller) are Cl (the default), Br, F, and I
        ion_concentration (float or u.Scalar[molarity]): ionic concentration in addition to
            whatever is needed to neutralize the solute. (if float is passed, we assume the
            number is Molar)

    Returns:
        moldesign.Molecule: new Molecule object containing both solvent and solute
    """
    import pdbfixer

    if padding is None and min_box_size is None:
        raise ValueError('Solvate arguments: must pass padding or min_box_size or both.')

    # add +s and -s to ion names if not already present
    if positive_ion[-1] != '+':
        assert positive_ion[-1] != '-'
        positive_ion += '+'
    if negative_ion[-1] != '-':
        assert negative_ion[-1] != '+'
        negative_ion += '-'

    if ion_concentration is not None:
        ion_concentration = u.MdtQuantity(ion_concentration)
        if ion_concentration.dimensionless:
            ion_concentration *= u.molar

    # calculate box size - in each dimension, use the largest of min_box_size or
    #    the calculated padding
    boxsize = np.zeros(3) * u.angstrom
    if min_box_size:
        boxsize[:] = min_box_size
    if padding:
        ranges = mol.positions.max(axis=0) - mol.positions.min(axis=0)
        for idim, r in enumerate(ranges):
            boxsize[idim] = max(boxsize[idim], r+padding)
    assert (boxsize >= 0.0).all()

    modeller = opm.mol_to_modeller(mol)

    # TODO: this is like 10 bad things at once. Should probably submit a
    #       PR to PDBFixer to make this a public staticmethod instead of a private instancemethod
    #       Alternatively, PR to create Fixers directly from Topology objs
    ff = pdbfixer.PDBFixer.__dict__['_createForceField'](None, modeller.getTopology(), True)

    modeller.addSolvent(ff,
                        boxSize=opm.pint2simtk(boxsize),
                        positiveIon=positive_ion,
                        negativeIon=negative_ion,
                        ionicStrength=opm.pint2simtk(ion_concentration),
                        neutralize=neutralize)

    return opm.topology_to_mol(modeller.getTopology(),
                               positions=modeller.getPositions(),
                               name='%s with water box' % mol.name)


