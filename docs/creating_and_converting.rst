Creating and converting molecules
---------------------------------


From other objects
==================
You can create a new molecule from any collection of atoms.

For instance, a list of atoms:
    >>> mol = mdt.Molecule([atom1, atom2])

An amino acid residue from another molecule:
   >>> protein = mdt.from_pdb('3AID')
   >>> mol = mdt.Molecule( protein.residues[3] )

Or even a list of molecules, atoms, and residues:
   >>> protein = mdt.from_pdb('3AID')
   >>> dmso = mdt.from_name('dmso')
   >>> cobalt_atom = mdt.Atom(symbol='Co')
   >>> complex = mdt.Molecule([protein, dmso, cobalt_atom])