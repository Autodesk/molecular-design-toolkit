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
""" PDB file parsing utilities

We don't yet (and hopefully will never need) an internal PDB parser or writer. For now,
routines in this module read and write data that's not necessarily parsed by other implementations.
"""
from past.builtins import basestring
import collections
from io import StringIO

import numpy as np

import moldesign as mdt

BioAssembly = collections.namedtuple('BioAssembly', 'desc chains transforms')


def get_conect_pairs(mol):
    """ Returns a dicitonary of HETATM bonds for a PDB CONECT record

    Note that this doesn't return the text records themselves, because they need
    to reference a specific PDB sequence number
    """
    conects = collections.OrderedDict()

    for residue in mol.residues:
        # intra-residue bonds
        if not residue.is_standard_residue:
            for bond in residue.bonds:
                if bond.order <= 1:
                    order = 1
                else:
                    order = bond.order
                for i in range(order):
                    conects.setdefault(bond.a1, []).append(bond.a2)

        # inter-residue bonds
        try:
            r2 = residue.next_residue
        except (StopIteration, KeyError, NotImplementedError):
            continue
        if not (residue.is_standard_residue and r2.is_standard_residue):
            for bond in residue.bonds_to(r2):
                conects.setdefault(bond.a1, []).append(bond.a2)

    return conects


def warn_assemblies(mol, assemblies):
    """ Print a warning message if the PDB structure contains a biomolecular assembly
    """

    # Don't warn if the only assembly is the asymmetric unit
    if len(assemblies) > 1 or len(list(assemblies.values())[0].transforms) > 1:
        print("WARNING: This PDB file contains the following biomolecular assemblies:")

        for name, asm in assemblies.items():
            print('WARNING: Assembly "%s": %d copies of chains %s'%(
                name, len(asm.transforms), ', '.join(asm.chains)))
            print('WARNING: Use ``mdt.build_assembly([molecule],[assembly_name])``' \
                  ' to build one of the above assemblies')


class MissingResidue(object):
    type = 'protein'
    missing = True

    def __init__(self, chain, resname, pdbindex):
        self.chain = chain
        self.resname = resname
        self.pdbindex = pdbindex

    @property
    def code(self):
        """str: one-letter amino acid code or two letter nucleic acid code, or '?' otherwise"""
        return mdt.data.RESIDUE_ONE_LETTER.get(self.resname, '?')


def get_pdb_missing_residues(fileobj):
    """ Parses missing residues from a PDB file.

    Args:
        fileobj (filelike): file-like access to the PDB file

    Returns:
        dict: listing of the missing residues of the form
           ``{[chain_id]:{[residue_number]:[residue_name], ...}, ...}``

    Examples:
        >>> with open('2jaj.pdb') as pdbfile:
        >>>     res = get_pdb_missing_residues(pdbfile)
        >>> res
        'A': {-4: 'GLY',
              -3: 'PRO',
              [...]
              284: 'SER'},
        'B': {34: 'GLY',
              35: 'GLU',
              [...]
    """
    lineiter = iter(fileobj)
    while True:
        try:
            fields = next(lineiter).split()
        except StopIteration:
            missing = {}  # no missing residues found in file
            break

        if fields == ['REMARK', '465', 'M', 'RES', 'C', 'SSSEQI']:
            missing = _parse_missing_xtal(fileobj)
            break
        elif fields[:3] == ['REMARK', '465', 'MODELS']:
            missing = _parse_missing_nmr(fileobj)
            break

    summary = {}
    for m in missing:
        summary.setdefault(m.chain, {})[m.pdbindex] = m.resname
    return summary


def _parse_missing_xtal(fileobj):
    missing = []
    while True:
        fields = next(fileobj).split()
        if fields[:2] != ['REMARK', '465']:
            break

        if len(fields) == 6:
            has_modelnum = 1
            if fields[2] != 1:  # only process the first model
                continue
        else:
            has_modelnum = 0

        missing.append(MissingResidue(chain=fields[3+has_modelnum],
                                      resname=fields[2+has_modelnum],
                                      pdbindex=int(fields[4+has_modelnum])))
    return missing


def _parse_missing_nmr(fileobj):
    header = next(fileobj).split()
    assert header == ['REMARK', '465', 'RES', 'C', 'SSSEQI']

    missing = []
    while True:
        fields = next(fileobj).split()

        if fields[:2] != ['REMARK', '465']:
            break

        missing.append(MissingResidue(chain=fields[3],
                                      resname=fields[2],
                                      pdbindex=int(fields[4])))
    return missing


def get_pdb_assemblies(fileobj):
    """Parse a PDB file, return biomolecular assembly specifications

    Args:
        fileobj (file-like): File-like object for the PDB file
            (this object will be rewound before returning)

    Returns:
        Mapping[str, BioAssembly]: dict mapping assembly ids to BioAssembly instances
    """
    assemblies = {}
    lineiter = iter(fileobj)
    while True:  # first, search for assembly transformations
        line = next(lineiter)
        fields = line.split()

        # Conditions that indicate we're past the "REMARK 350" section
        if fields[0] in ('ATOM', 'HETATM', 'CONECT'):
            break
        if fields[0] == 'REMARK' and int(fields[1]) > 350:
            break

        # look for start of a assembly transformation, i.e. "REMARK 350 BIOMOLECULE: [name]  "
        if fields[:3] == 'REMARK 350 BIOMOLECULE:'.split():
            assembly_name = fields[-1]
            assemblies[assembly_name] = _read_pdb_assembly(lineiter)

    return assemblies


def _read_pdb_assembly(lineiter):
    """Helper for get_pdb_assemblies
    """
    # First, there's description lines: "REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: OCTAMERIC"
    description_lines = []
    line = next(lineiter)
    fields = line.split()
    while fields[:7] != 'REMARK 350 APPLY THE FOLLOWING TO CHAINS:'.split():
        description_lines.append(line[len('REMARK 350 '):])
        line = next(lineiter)
        fields = line.split()
    description = (''.join(description_lines)).strip()

    # Next, we get the chains in this assembly: "REMARK 350 APPLY THE FOLLOWING TO CHAINS:  C, D"
    assert fields[:7] == 'REMARK 350 APPLY THE FOLLOWING TO CHAINS:'.split()
    chain_names = [x.rstrip(',') for x in fields[7:]]
    while fields[-1][-1] == ',':  # deal with multi-line lists of chains
        line = next(lineiter)
        fields = line.split()
        assert fields[2:4] == ['AND', 'CHAINS:']
        chain_names.extend(x.rstrip(',') for x in fields[4:])

    transforms = []
    while True:  # loop over each assembly transformation
        # example: "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000 "
        t = np.zeros((4, 4))
        t[3, 3] = 1.0

        for idim in range(3):
            line = next(lineiter)
            fields = line.split()
            if idim == 0 and len(fields) == 2:
                return BioAssembly(description, chain_names, transforms)

            assert int(fields[3]) == len(transforms)+1
            assert fields[2] == ('BIOMT%d' % (idim+1))
            t[idim, :] = list(map(float, fields[4:8]))

        transforms.append(t)


def assign_biopolymer_bonds(mol):
    """ Assign bonds to all standard residues using the PDB chemical component dictionary

    Any unrecognized residues are ignored.

    References:
        http://www.wwpdb.org/data/ccd
    """
    for chain in mol.chains:
        try:
            chain.assign_biopolymer_bonds()
        except KeyError:
            print(('WARNING: failed to assign backbone bonds for %s') % str(chain))
    for residue in mol.residues:
        try:
            residue.assign_template_bonds()
        except KeyError:
            if residue.type not in ('ion', 'water'):
                print(('WARNING: failed to assign bonds for %s; use '
                      '``residue.assign_distance.bonds`` to guess the topology') % str(residue))


def assign_unique_hydrogen_names(mol):
    """ Assign unique names to all hydrogens, based on either:
    1) information in the Chemical Component Database, or
    2) newly generated, unique names

    Args:
        mol (moldesign.Molecule):
    """

    for residue in mol.residues:
        if residue.resname in mdt.data.RESIDUE_BONDS:
            _assign_hydrogen_names_from_ccd(residue)
        else:
            _assign_unique_hydrogen_names_in_order(residue)

        residue.rebuild()


def _assign_hydrogen_names_from_ccd(residue):
    ccd_bonds = mdt.data.RESIDUE_BONDS[residue.resname]
    taken = set(atom.name for atom in residue.atoms)

    if 'H' not in taken:
        return  # nothing to do
    if 'H' in ccd_bonds:
        taken.remove('H')  # someone will actually need to be named "H'

    for atom in residue:
        if atom.atnum != 1 or atom.name != 'H':
            continue
        assert atom.num_bonds == 1, 'Hydrogen has more than one bond'
        bond = atom.bonds[0]
        other = bond.partner(atom).name
        for partner in ccd_bonds[other]:
            if partner[0] == 'H' and partner not in taken:
                assert ccd_bonds[other][partner] == 1, 'Hydrogen bond order is not 1'
                atom.name = partner
                taken.add(partner)
                break


def _assign_unique_hydrogen_names_in_order(residue):
    n_hydrogen = 1
    namecounts = collections.Counter(x.name for x in residue.atoms)
    if namecounts.get('H', 0) > 1:
        used_names = set(atom.name for atom in residue.atoms)
        for atom in residue.atoms:
            if atom.name == 'H':
                name = 'H%d' % n_hydrogen
                while name in used_names:
                    n_hydrogen += 1
                    name = 'H%d' % n_hydrogen
                atom.name = name
                used_names.add(name)
