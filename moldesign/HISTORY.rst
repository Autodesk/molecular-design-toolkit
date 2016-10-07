0.7.3 (WIP)
===========
NEW MODELING FEATURES
 - #33 - Add DFT w/ gradients; MP2, CASSCF, CASCI w/out gradients
 - Constrained minimizations w/ SHAKE and scipy's SLQSP
 - Rename `build_bdna` to `build_dna_helix` and give access to all NAB helix types
 - GAFF parameterizer for small molecules -- ``params = mdt.parameterize(mol)``
 - AM1-BCC and Gasteiger partial charge calculators: ``mdt.calc_am1_bcc_charges`` and
    ``mdt.calc_gasteiger_charges``
 - Add PDB database and biomolecular assembly support for mmCIF files
 - #72 - Add ``moldesign.guess_formal_charges`` and ``moldesign.add_missing_data``
 - Excited and multi-state property calculations

OTHER ENHANCEMENTS
 - #78 - ``moldesign`` now imports much more quickly
 - Add ``GAFF`` energy model automates small molecule parameterization
 - Change Example 2 to show an absorption spectrum calculation
 - Add Example 4 on protein MD with a small ligand
 - Add Example 5: on constrained minimization and enthalpic barriers
 - Add Tutorial 3: QM data analysis
 - Show changelog and version check in the ``mdt.about()`` (aka ``mdt.configure``) widget
 - Change moldesign.tools and moldesign.helpers modules into more rationally organized subpackages
 - ``mdt.set_dihedral`` can be called with two atoms in the same way as ``mdt.dihedral``
 - Explicit parameter created to store wavefunction guesses
 - Better access to density matrix in wavefunction objects
 - Improved parsing support for PDB and mmCIF files

BUGFIXES
 - #61 - fixed a KeyError when parsing PDBs with metal centers or ions
 - #74 - Add function to create PDB files with correct TER records (used for TLeap input)
 - Better handling of chains with non-standard residues
 - ``mdt.add_hydrogens`` no longer creates structures with separated residues
 - Fix sign of dihedral gradient
 - Charge quantities now mostly have the correct units


0.7.2 - July 26, 2016
=====================
NEW MODELING FEATURES
 - Add trajectory geometry analysis functions (``traj.dihedral``, ``traj.angle``, etc.)
 - Can now calculate angles, dists, and dihedrals by name within a residue
    (``residue.distance('CA','CG')``)
 - Calculate dihedral angles using only two atoms defining the central bond (MDT will infer
    infer the other two in a consistent way)

CHANGES
 - Completed tutorials

BUGFIXES
 - #28: Fixed a rounding error and logging problems with OpenMM trajectory snapshots
 - #21: Better bond orders to structures in Amber files, which don't store them
 - #20: Store OpenMM force vector correctly

0.7.1 - July 20, 2016
=====================
BUGFIXES
  - #4: Use public demo CCC server by default
  - #3: Fix ``python -m moldesign intro``

0.7.0 - July 15, 2016
=====================
 - Initial public release
