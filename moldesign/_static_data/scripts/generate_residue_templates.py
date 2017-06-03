#!/usr/bin/env python
""" A short script that assembles reasonable template geometries for standard amino acids
and DNA bases
"""
import moldesign as mdt
import json

protein = mdt.from_pdb('4F0A')
dna = mdt.build_bdna('ACTGAA')

residues = {}

for chain in list(protein.chains.values()) + list(dna.chains.values()):
    for residue in chain:
        if residue.type not in ('dna', 'protein'): continue
        if residue.type == 'dna' and (residue.is_5prime_end or residue.is_3prime_end):
            continue
        if residue.type == 'protein' and (residue.is_c_terminal or residue.is_n_terminal):
            continue

        residues[residue.resname] = mdt.clean_pdb(mdt.Molecule(residue))

residue_pdb = {k: m.write(format='pdb') for k,m in residues.items()}

with open('residue_templates.json', 'w') as outfile:
    json.dump(residue_pdb, outfile)
