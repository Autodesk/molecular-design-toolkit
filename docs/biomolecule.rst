Biomolecular structure
======================



Primary structure
-----------------
**Class documentation:** :class:`moldesign.Chain`, :class:`moldesign.Residue`

Biomolecules also contain primary structure information such as :class:`Chains <moldesign.Chain>`
and :class:`Residues <moldesign.Residue>`. Chains can be accessed by name OR by index:

   >>> chain1 = molecule.chains['A']
   >>> chain2 = molecule.chains[0]
   >>> chain1 is chain2
   True

Each chain contains :class:`residues <moldesign.Residue>`. In a chain, residues can similarly be
accessed through a flat list or by name:

   >>> res0 = molecule.residues[0]
   >>> resA = molecule.chains['A'].residues['PRO1']
   >>> res0 is resA
   True

A flat list of all residues in a molecule is also available at `molecule.residues`.



Biomolecular assemblies
-----------------------
Many biomolecules in the PDB only contain a subset of the total biomolecular structure - the
remaining parts of the structure can be generated via `symmetry transformations <http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies>`_.

When you read in such a structure, MDT will issue a warning.

   >>> mol = mdt.from_pdb('3FPP')
   WARNING: This PDB file contains the following biomolecular assemblies:
   WARNING: Assembly "1": 3 copies of chains A, B
   WARNING: Use ``mdt.build_assembly([molecule],[assembly_name])`` to build one of the above assemblies

To create the full assembly, run
   >>> assembly = mdt.build_assembly(mol,"1")
   >>> assembly.draw()

   .. image:: img/howdoi_pdb_assm.png

Note:
   Only PDB-formatted files are currently supported for biomolecular assemblies - MMCif support
   is in progress.





