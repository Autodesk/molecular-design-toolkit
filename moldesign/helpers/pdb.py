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
from cStringIO import StringIO

import numpy as np

import moldesign as mdt

BioAssembly = collections.namedtuple('BioAssembly', 'desc chains transforms')


def insert_ter_records(mol, pdbfile):
    """ Inserts TER records to indicate the end of the biopolymeric part of a chain

    Many common PDB writers - including OpenBabel - don't insert TER records. This can
    cause a problem for situations such as forcefield assignment. This routine is one
    solution to that problem.

    What it does:
    This routine inserts 'TER' records
    1) after any protein residue not bound to the next residue via backbone peptide bond, and
    2) after any DNA residue not bound to the next residue via a backbone phosphate bond

    In the input PDB file, the ATOM records be formatted with the proper columns (see
    http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM) and
    the chain names and residue numbers must match ``chain.pdbname`` and ``residue.pdbindex``,
    respectively.

    Args:
        mol (moldesign.Molecule): the MDT version of the molecule that pdbfile describes
        pdbfile (TextIO OR str): pdb file (file-like or string)

    Returns:
        cStringIO.StringIO OR str: copy of the pdb file with added TER records - it will be
         returned as the same type passed (i.e., as a filelike buffer or as a string)
    """

    is_string = False
    if isinstance(pdbfile, basestring):
        pdbfile = StringIO(pdbfile)
        is_string = True

    # First, identify where to insert records (if anywhere)
    ter_residues = set()
    for chain in mol.chains:
        if chain.type == 'protein':
            ter_residues.add((chain.pdbname, chain.c_terminal.pdbindex))
        elif chain.type == 'dna':
            ter_residues.add((chain.pdbname, chain.threeprime_end.pdbindex))

    # insert the records (if necessary)
    newf = StringIO()
    pdbfile.seek(0)
    watchres = None
    for line in pdbfile:
        fields = line.split()
        if fields and fields[0] in ('ATOM','HETATM'):  # line is an atom record
            res = (line[21], int(line[22:26].strip()))
            if watchres:
                if res != watchres:
                    print >> newf, 'TER'
                    watchres = None
            if res in ter_residues:
                watchres = res

        elif watchres is not None:  # line is not an atom record
            watchres = None
            if line.strip() != 'TER':
                print >> newf, 'TER'

        newf.write(line)

    newf.seek(0)
    if is_string:
        return newf.read()
    else:
        return newf

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


def guess_atnum_from_name(s):
    """ Guess an atomic number given a name string (usually 1-3 characters).

    Args:
        s (str): atomic number

    Returns:
        int: atomic number

    Raises:
        KeyError: if atomic number can't be determined

    Examples:
        >>> guess_atnum_from_name('C')
        6
        >>> guess_atnum_from_name('C1')
        6
        >>> guess_atnum_from_name('cl3')
        17
        >>> guess_atnum_from_name('CL')
        17
    """
    try:  # the unmodified string
        return mdt.data.ATOMIC_NUMBERS[s]
    except KeyError:
        pass

    cleaned = ''.join((c.upper() if i==0 else c.lower())
                      for i,c in enumerate(s)
                      if c.isalpha())

    try:  # just the letters, with proper capitalization
        return mdt.data.ATOMIC_NUMBERS[cleaned]
    except KeyError:
        pass

    # otherwise, just the first letter
    return mdt.data.ATOMIC_NUMBERS[cleaned[0]]


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
        if len(fields) == 0:
            continue
        if fields[0] != 'CONECT':
            continue

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
        try:
            chain.assign_biopolymer_bonds()
        except KeyError:
            print('WARNING: failed to assign backbone bonds for %s') % str(chain)
    for residue in mol.residues:
        try:
            residue.assign_template_bonds()
        except KeyError:
            print('WARNING: failed to assign bonds for %s; use '
                  '``residue.assign_distance.bonds`` to guess the topology') % str(residue)


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
