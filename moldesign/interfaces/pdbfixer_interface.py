from __future__ import print_function, absolute_import, division

import string

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
from past.builtins import basestring

import imp
import io
import re


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
            pdbfile=io.StringIO(mol.write(format='pdb')))
    return fixer


@exports
def fixer_to_mol(f):
    return opm.topology_to_mol(f.topology, positions=f.positions)


@compute.runsremotely(enable=force_remote)
def mutate_residues(mol, residue_map):
    """ Create a mutant with point mutations (returns a copy - leaves the original unchanged)

    Mutations may be specified in one of two ways:
      1) As a dictionary mapping residue objects to the 3-letter name of the amino acid they
        will be mutated into: ``{mol.residues[3]: 'ALA'}``
      2) As a list of mutation strings: ``['A43M', '332S', 'B.53N']`` (see below)

    Mutation strings have the form:
       ``[chain name .][initial amino acid code](residue index)(mutant amino acid code)``
       If the chain name is omited, the mutation will be applied to all chains (if possible).
       The initial amino acid code may also be omitted.

    Examples:
        >>> mutate_residues(mol, {mol.residues[5]: 'ALA'})  # mutate residue 5 to ALA
        >>> mutate_residues(mol, 'A43M')  # In all chains with ALA43 mutate it to MET43
        >>> mutate_residues(mol, ['A.332S', 'B.120S'])  # Mutate Chain A res 332 and B 120 to SER
        >>> mutate_residues(mol, ['B.C53N']) # Mutate Cysteine 53 in chain B to Asparagine

    Args:
        mol (moldesign.Molecule): molecule to create mutant from
        residue_map (Dict[moldesign.Residue:str] OR List[str]): list of mutations (see above for
           allowed formats)

    Returns:
        moldesign.Molecule: the mutant
    """
    fixer = mol_to_fixer(mol)
    chain_mutations = {}
    mutation_strs = []
    if not hasattr(residue_map, 'items'):
        residue_map = _mut_strs_to_map(mol, residue_map)

    if not residue_map:
        raise ValueError("No mutations specified!")

    for res, newres in residue_map.items():
        chain_mutations.setdefault(res.chain.pdbname, {})[res] = residue_map[res]
        mutation_strs.append(_mutation_as_str(res, newres))

    for chainid, mutations in chain_mutations.items():
        mutstrings = ['%s-%d-%s' % (res.resname, res.pdbindex, newname)
                      for res, newname in mutations.items()]
        fixer.applyMutations(mutstrings, chainid)
    temp_mutant = fixer_to_mol(fixer)
    _pdbfixer_chainnames_to_letter(temp_mutant)

    # PDBFixer reorders atoms, so to keep things consistent, we'll graft the mutated residues
    # into an MDT structure
    assert temp_mutant.num_residues == mol.num_residues  # shouldn't change number of residues
    residues_to_copy = []
    for oldres, mutant_res in zip(mol.residues, temp_mutant.residues):
        if oldres in residue_map:
            residues_to_copy.append(mutant_res)
            mutant_res.mol = None
            mutant_res.chain = oldres.chain
            for atom in mutant_res.atoms:
                atom.chain = oldres.chain
        else:
            residues_to_copy.append(oldres)

    metadata = {'origin': mol.metadata.copy(),
                'mutations': mutation_strs}

    return mdt.Molecule(residues_to_copy,
                        name='Mutant of "%s"' % mol,
                        metadata=metadata)


def _pdbfixer_chainnames_to_letter(pdbfixermol):
    for chain in pdbfixermol.chains:
        try:
            if chain.name.isdigit():
                chain.name = string.ascii_uppercase[int(chain.name)-1]
        except (ValueError, TypeError, IndexError):
            continue  # not worth crashing over


def _mutation_as_str(res, newres):
    """ Create mutation string for storage as metadata.

    Note that this will include the name of the chain, if available, as a prefix:

    Examples:
        >>> res = mdt.Residue(resname='ALA', pdbindex='23', chain=mdt.Chain(name=None))
        >>> _mutation_as_str(res, 'TRP')
        'A23W'
        >>> res = mdt.Residue(resname='ALA', pdbindex='23', chain=mdt.Chain(name='C'))
        >>> _mutation_as_str(res, 'TRP')
        'C.A23W'

    Args:
        res (moldesign.Residue): residue to be mutated
        newres (str): 3-letter residue code for new amino acid

    Returns:
        str: mutation string

    References:
        Nomenclature for the description of sequence variations
            J.T. den Dunnen, S.E. Antonarakis: Hum Genet 109(1): 121-124, 2001
            Online at http://www.hgmd.cf.ac.uk/docs/mut_nom.html#protein
    """
    try:  # tries to describe mutation using standard
        mutstr = '%s%s%s' % (res.code, res.pdbindex,
                             mdt.data.RESIDUE_ONE_LETTER.get(newres, '?'))
        if res.chain.pdbname:
            mutstr = '%s.%s' % (res.chain.pdbname, mutstr)
        return mutstr
    except (TypeError, ValueError) as e:
        print('WARNING: failed to construct mutation code: %s' % e)
        return '%s -> %s' % (str(res), newres)


MUT_RE = re.compile(r'(.*\.)?([^\d]*)(\d+)([^\d]+)')  # parses mutation strings


def _mut_strs_to_map(mol, strs):
    if isinstance(strs, basestring):
        strs = [strs]
    mutmap = {}
    for s in strs:
        match = MUT_RE.match(s)
        if match is None:
            raise ValueError("Failed to parse mutation string '%s'" % s)
        chainid, initialcode, residx, finalcode = match.groups()
        if chainid is not None:
            parent = mol.chains[chainid[:-1]]
        else:
            parent = mol  # queries the whole molecule

        newresname = mdt.data.RESIDUE_CODE_TO_NAME[finalcode]

        query = {'pdbindex': int(residx)}
        if initialcode:
            query['code'] = initialcode

        residues = parent.get_residues(**query)

        if len(residues) == 0:
            raise ValueError("Mutation '%s' did not match any residues" % s)

        for res in residues:
            assert res not in mutmap, "Multiple mutations for %s" % res
            mutmap[res] = newresname

    return mutmap


@compute.runsremotely(enable=force_remote)
def add_water(mol, min_box_size=None, padding=None,
              ion_concentration=0.0, neutralize=True,
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

    ion_concentration = u.MdtQuantity(ion_concentration)
    if ion_concentration.dimensionless:
        ion_concentration *= u.molar
    ion_concentration = opm.pint2simtk(ion_concentration)

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

    # Creating my fixers directly from Topology objs
    ff = pdbfixer.PDBFixer.__dict__['_createForceField'](None, modeller.getTopology(), True)

    modeller.addSolvent(ff,
                        boxSize=opm.pint2simtk(boxsize),
                        positiveIon=positive_ion,
                        negativeIon=negative_ion,
                        ionicStrength=ion_concentration,
                        neutralize=neutralize)

    solv_tempmol = opm.topology_to_mol(modeller.getTopology(),
                                  positions=modeller.getPositions(),
                                  name='%s with water box' % mol.name)

    _pdbfixer_chainnames_to_letter(solv_tempmol)


    # PDBFixer reorders atoms, so to keep things consistent, we'll graft the mutated residues
    # into an MDT structure
    newmol_atoms = [mol]
    for residue in solv_tempmol.residues[mol.num_residues:]:
        newmol_atoms.append(residue)

    newmol = mdt.Molecule(newmol_atoms,
                          name="Solvated %s" % mol,
                          metadata={'origin':mol.metadata})

    assert newmol.num_atoms == solv_tempmol.num_atoms
    assert newmol.num_residues == solv_tempmol.num_residues
    return newmol



