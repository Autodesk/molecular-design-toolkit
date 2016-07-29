0.7.3 (WIP)
===========
NEW MODELLING FEATURES
 - GAFF parameterizer for small molecules -- ``params = mdt.parameterize(mol)``
 - AM1-BCC and Gasteiger partial charge calculators

CHANGES
 - create changelog and version check to the ``mdt.about()`` (aka mdt.configure) widget

BUGFIXES


0.7.2 - July 26, 2016
=====================
NEW MODELLING FEATURES
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
