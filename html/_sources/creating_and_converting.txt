Reading and building molecules
==============================


From names and IDs
------------------
You can use an IUPAC name:

   >>> import moldesign as mdt
   >>> benzene = mdt.from_name('benzene')
   >>> caffeine = mdt.from_name('1,3,7-Trimethylpurine-2,6-dione')

or a `SMILES string <https://en.wikipedia.org/wiki/Simplified_molecular-input_line-entry_system>`_

   >>> benzene = mdt.from_smiles('c1ccccc1')
   >>> caffeine = mdt.from_smiles('CN1C=NC2=C1C(=O)N(C(=O)N2C)C')


To download a biomolecular structure from the `Protein DataBank <http://www.rcsb.org/>`_, you can use

   >>> hiv_protease = mdt.from_pdb('3AID')

DNA helices can be generated from a genetic sequence:

   >>> bdna = mdt.build_dna_helix('ACTG')


From files
----------
MDT supports most common molecular formats via ``moldesign.read``. In most cases, the format can be inferred from the filename. Files compressed with `gz` or `bz2` can be read as well.

   >>> benzene = mdt.read('benzene.sdf')
   >>> caffeine = mdt.read('caffeine.xyz')
   >>> lsd = mdt.read('lsd.mol2.gz')
   >>> hiv_protease = mdt.read('3aid.pdb.bz2')
   >>> zika_capsid = mdt.read('5ire.cif')

In addition, any pickled python object can be read in with this function - files with ``pickle``, ``pkl``, and ``P`` are recognized as pickle files:
   >>> saved_mol = mdt.read('molecule_1.pickle')
   >>> saved_atom = mdt.read('atom_2.pkl.gz')
   >>> saved_trajectory = mdt.read('traj.P.bz2')


From strings
------------
File content in strings or file-like objects can be read as well, but the format needs to be explicitly specified.

   >>> water_str = """ 3
   >>>                 water xyz file
   >>>                 O          0.98285        0.07497        0.04837
   >>>                 H          0.70400        0.94631        0.36769
   >>>                 H          1.95074        0.11856        0.06434 """
   >>> water = mdt.read(water_str, format='xyz')
   >>>
   >>> import StringIO
   >>> water_filelike = StringIO.StringIO(water_str)
   >>> molecule = mdt.read(water_filelike, format='xyz')


From other molecules
--------------------
You can create a new molecule from any collection of atoms.

For instance, a list of atoms:
    >>> mol = mdt.Molecule([atom1, atom2])

An amino acid residue from another molecule:
   >>> protein = mdt.from_pdb('3AID')
   >>> mol1 = mdt.Molecule(protein.atoms[0:20])
   >>> mol2 = mdt.Molecule(protein.chains['A'].residue['PRO1'])

Or even a list of molecules, atoms, and residues:
   >>> protein = mdt.from_pdb('3AID')
   >>> dmso = mdt.from_name('dmso')
   >>> cobalt_atom = mdt.Atom(symbol='Co')
   >>> complex = mdt.Molecule([protein, dmso, cobalt_atom])


From other python packages
--------------------------
MDT's interfaces allow it to import objects from a variety of other molecular modeling packages, including;

- `OpenBabel <http://openbabel.org/>`_ / `Pybel <https://openbabel.org/docs/dev/UseTheLibrary/Python_Pybel.html>`_ : ``mdt.interfaces.pybel_to_mdt(pybel_molecule)``
- `OpenMM <http://openmm.org/>`_: ``mdt.interfaces.openmm_topology_to_mdt(openmm_topology)``
- `BioPython <http://biopython.org/wiki/Biopython>`_: ``mdt.interfaces.biopython_to_mdt(biopython_pdb_structure)``
- `PdbFixer <https://github.com/pandegroup/pdbfixer>`_: ``mdt.interfaces.pdbfixer_to_mdt(pdbfixer_object)``