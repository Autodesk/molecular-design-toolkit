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
""" PDB file parsing utilities

We don't yet (and hopefully will never need) an internal PDB parser or writer. For now,
routines in this module read and write data that's not necessarily parsed by other implementations.
"""
import collections

import numpy as np

BioAssembly = collections.namedtuple('BioAssembly', 'desc chains transforms')


def warn_assemblies(mol, assemblies):
    """ Print a warning message if the PDB structure contains a biomolecular assembly
    """

    # Don't warn if the only assembly is the asymmetric unit
    if len(assemblies) > 1 or len(assemblies.values()[0].transforms) > 1:
        print "WARNING: This PDB file contains the following biomolecular assemblies:"
        for name, asm in assemblies.iteritems():
            print 'WARNING: Assembly "%s": %d copies of chains %s'%(
                name, len(asm.transforms), ', '.join(asm.chains))
            print 'WARNING: Use ``mdt.build_assembly([molecule],[assembly_name])``' \
                  ' to build one of the above assemblies'


def get_conect_records(pdbfile):
    """Parse a PDB file, return CONECT records

    Bond orders are assigned as 1 by default. Repeated CONECT records are interpreted as
    higher order bonds.

    Args:
        pdbfile (file): file-like object

    Example:
        > CONECT   1   2   3
        > CONECT   1   2
        > CONECT   2   1   1
        > CONECT   3   1
        These records are interpreted as a double-bond between atoms 1 and 2
        and a single bond between atoms 1 and 3

    Note:
        This only returns the covalent CONECT records (the first 4 entries) - it doesn't
        return salt bridges or hydrogen bonds

    Returns:
        dict: CONECT records using serial numbers -
             ``{serial1: {serial2:order. serial3:order, }, ...}``
    """
    conect = {}
    for line in pdbfile:
        fields = line.split()
        if len(fields) == 0: continue
        if fields[0] != 'CONECT': continue

        atombonds = conect.setdefault(int(fields[1]), {})
        for f in fields[2:6]:  # TODO: check the end bound
            serial = int(f)
            if serial not in atombonds:
                atombonds[serial] = 0
            atombonds[serial] += 1
    return conect


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
        line = lineiter.next()
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
    line = lineiter.next()
    fields = line.split()
    while fields[:7] != 'REMARK 350 APPLY THE FOLLOWING TO CHAINS:'.split():
        description_lines.append(line[len('REMARK 350 '):])
        line = lineiter.next()
        fields = line.split()
    description = (''.join(description_lines)).strip()

    # Next, we get the chains in this assembly: "REMARK 350 APPLY THE FOLLOWING TO CHAINS:  C, D"
    assert fields[:7] == 'REMARK 350 APPLY THE FOLLOWING TO CHAINS:'.split()
    chain_names = [x.rstrip(',') for x in fields[7:]]
    while fields[-1][-1] == ',':  # deal with multi-line lists of chains
        line = lineiter.next()
        fields = line.split()
        assert fields[2:4] == ['AND', 'CHAINS:']
        chain_names.extend(x.rstrip(',') for x in fields[4:])

    transforms = []
    while True:  # loop over each assembly transformation
        # example: "REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000 "
        t = np.zeros((4, 4))
        t[3, 3] = 1.0

        for idim in xrange(3):
            line = lineiter.next()
            fields = line.split()
            if idim == 0 and len(fields) == 2:
                return BioAssembly(description, chain_names, transforms)

            assert int(fields[3]) == len(transforms)+1
            assert fields[2] == ('BIOMT%d' % (idim+1))
            t[idim, :] = map(float, fields[4:8])

        transforms.append(t)


def assign_biopolymer_bonds(mol):
    """ Assign bonds to all standard residues using the PDB chemical component dictionary

    Any unrecognized residues are ignored.

    References:
        http://www.wwpdb.org/data/ccd
    """

    for chain in mol.chains:
        chain.assign_biopolymer_bonds()
    for residue in mol.residues:
        try:
            residue.assign_template_bonds()
        except KeyError:
            print ('WARNING: failed to assign bonds for residue %s; use '
                   '``residue.assign_distance.bonds`` to guess the topology') % str(residue)