0.7.3 (WIP)
===========
 - Add changelog and version check to the ``mdt.about()`` (aka mdt.configure) widget
 - Add GAFF parameterizer for small molecules -- ``params = mdt.parameterize(mol)``


0.7.2 - July 26, 2016
=====================
FEATURES
 - Completed tutorials
 - Add trajectory geometry analysis functions (``traj.dihedral``, ``traj.angle``, etc.)
 - Can now calculate angles, dists, and dihedrals by name within a residue
    (``residue.distance('CA','CG')``)
 - Can request a dihedral angle using only two atoms defining the central bond (MDT will try to
    infer the other two atoms uniquely and sanely)

BUGS
 - #28: Fixed a rounding error and logging problems with OpenMM trajectory snapshots
 - #21: Better bond orders to structures in Amber files, which don't store them
 - #20: Store OpenMM force vector correctly

0.7.1 - July 20, 2016
=====================
Bugfix release

BUGS
  - #4: Use public demo CCC server by default
  - #3: Fix ``python -m moldesign intro``

0.7.0 - July 15, 2016
=====================
 - Initial public release
