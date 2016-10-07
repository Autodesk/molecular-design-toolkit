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

import string

import moldesign as mdt
import moldesign.molecules

from moldesign.interfaces.ambertools import build_bdna, build_dna_helix

from . import toplevel, __all__ as _pkgall

_pkgall.extend(['build_bdna', 'build_dna_helix'])


@toplevel
def build_assembly(mol, assembly_name):
    """ Create biological assembly using a bioassembly specification.

    This routine builds a biomolecular assembly using the specification from a PDB header (if
    present, this data can be found in the  "REMARK 350" lines in the PDB file). Assemblies are
    author-assigned structures created by copying, translating, and rotating a subset of the
    chains in the PDB file.

    See Also:
        http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies

    Args:
        mol (moldesign.Molecule): Molecule with assembly data (assembly data will be created by the
            PDB parser at ``molecule.properties.bioassembly``)
        assembly_name (str OR int): id of the biomolecular assembly to build.

    Returns:
        mol (moldesign.Molecule): molecule containing the complete assembly

    Raises:
        AttributeError: If the molecule does not contain any biomolecular assembly data
        KeyError: If the specified assembly is not present
    """
    if isinstance(assembly_name, int): assembly_name = str(assembly_name)

    if 'bioassemblies' not in mol.properties:
        raise AttributeError('This molecule does not contain any biomolecular assembly data')
    try:
        asm = mol.properties.bioassemblies[assembly_name]
    except KeyError:
        raise KeyError(('The specified assembly name ("%s") was not found. The following '
                        'assemblies are present: %s') %
                       (assembly_name,
                        ', '.join(mol.properties.bioassemblies.keys())))

    # Make sure each chain gets a unique name - up to all the letters in the alphabet, anyway
    used_chain_names = set()
    alpha = iter(string.ascii_uppercase)

    # Create the new molecule by copying, transforming, and renaming the original chains
    all_atoms = moldesign.molecules.atomcollections.AtomList()
    for i, t in enumerate(asm.transforms):
        for chain_name in asm.chains:
            chain = mol.chains[chain_name].copy()
            chain.transform(t)

            while chain.name in used_chain_names:
                chain.name = alpha.next()
            used_chain_names.add(chain.name)
            chain.pdbname = chain.pdbindex = chain.name
            all_atoms.extend(chain.atoms)
    newmol = mdt.Molecule(all_atoms,
                          name="%s (bioassembly %s)" % (mol.name, assembly_name))
    return newmol