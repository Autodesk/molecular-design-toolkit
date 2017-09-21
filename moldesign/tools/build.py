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
                        ', '.join(list(mol.properties.bioassemblies.keys()))))

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
                chain.name = next(alpha)
            used_chain_names.add(chain.name)
            chain.pdbname = chain.pdbindex = chain.name
            all_atoms.extend(chain.atoms)
    newmol = mdt.Molecule(all_atoms,
                          name="%s (bioassembly %s)" % (mol.name, assembly_name))
    return newmol

@toplevel
def build_model(temp_molecule, temp_chain_id, temp_start, temp_end, temp_sequence,
        targ_start, targ_end, targ_sequence): 
    """ Build a model of a structure for a given target sequence from a 
        template sequence and structure.

    Note that there is often a discrepancy between the number of residues in the temp
    sequence and the number of actual residues contained in the molecule. 

    Args:
        temp_molecule(moldesign.Molecule): The molecule used as a template structure. 
        temp_chain_id(string): The ID of the chain in the template structure. 
        temp_start(int): The start of the template sequence. 
        temp_end(int): The end of the template sequence. 
        temp_sequence(string): The 1 character amino acid sequence of the template.
        targ_start(int): The start of the target sequence. 
        targ_end(int): The end of the target sequence. 
        targ_sequence(string): The 1 charcter amino acid sequence of the target. 

    Returns:
        moldesign.Molecule: Molecule representing a structural model.
    """
    if temp_chain_id not in temp_molecule.chains:
        raise ValueError("No chain '%s' found in molecule '%s'" % 
             (temp_chain_id,temp_molecule.name))

    chain = temp_molecule.chains[temp_chain_id]
    residues = chain.residues
    atoms = chain.atoms

    # Add residues from the template sequence that are present in the template
    # molecule and create a list of mutations where they differ from the target 
    # sequence. 
    mutation_list = []
    temp_residues = []
    temp_seq_num = temp_start
    targ_seq_num = targ_start
    for temp_aa, targ_aa in zip(temp_sequence, targ_sequence):
        temp_res_name = mdt.data.mdt.data.RESIDUE_CODE_TO_NAME[temp_aa] + str(temp_seq_num)
        if temp_res_name in residues:
            if temp_aa != targ_aa:
                mutate_str = "%s.%s%d%s" % (temp_chain_id, temp_aa, targ_seq_num, targ_aa)
                mutation_list.append(mutate_str)
            temp_res = residues[temp_res_name]
            # Renumber the residue for the target sequence.
            temp_res.pdbindex = targ_seq_num 
            temp_res.name = mdt.data.mdt.data.RESIDUE_CODE_TO_NAME[temp_aa] + str(targ_seq_num)
            temp_residues.append(temp_res)
        temp_seq_num += 1
        targ_seq_num += 1
    #__for temp_aa, targ_aa in zip(temp_sequence, targ_sequence)

    # Create a template molecule from the template sequence. 
    temp_molecule = mdt.Molecule(temp_residues)
    if len(mutation_list) == 0:
        return temp_molecule

    # Change residues in template molecule to residues of the target sequence.
    model_molecule = mdt.mutate_residues(temp_molecule, mutation_list)

    # Check that the number of bonds for backbone atoms match.
    for res, model_res in zip(temp_molecule.residues, model_molecule.residues):
        if not res.backbone:
            continue
        for atom in res.backbone:
            bonds = [bond for bond in temp_molecule.bond_graph[atom] if bond.name in res.backbone]
            model_atom = model_molecule.chains[temp_chain_id].residues[model_res.name].atoms[atom.name]
            model_bonds = [bond for bond in model_molecule.bond_graph[model_atom] \
                if bond.name in model_res.backbone]
            if len(bonds) != len(model_bonds):
                raise Exception("Backbone atom bonds are missing")

    return model_molecule 

