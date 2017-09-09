## Changelog

### 0.8.0 - September 9, 2017

MDT 0.8 represents a substantial refinement of the core MDT code, offering Python 2/3 support
and increased stability and robustness.

##### NEW MODELING FEATURES
  - Initial **NWChem** integration and pre-compiled docker image ([\#120](https://github.com/autodesk/molecular-design-toolkit/issues/120))
  - **ParmEd integration** - assigned forcefield parameters are now stored as ParmEd objects, and
    MDT `Molecule` objects can be interconverted with ParmEd `Structures` via
    `mdt.interfaces.parmed_to_mdt` and `mdt.interfaces.mdt_to_parmed` ([\#116](https://github.com/autodesk/molecular-design-toolkit/issues/116))
  - **Orbital and basis function** descriptions are now mathematically complete, allowing
    operations such as gaussian multiplication, overlap calculations, and real-space amplitude
    evaluation ([\#167](https://github.com/autodesk/molecular-design-toolkit/issues/167))
  - Overhauled **forcefield handling** ([\#149](https://github.com/autodesk/molecular-design-toolkit/issues/149)): 
    1. `mdt.assign_forcefield` has been replaced with flexible `ForceField` objects, which offer
        `Forcefield.assign` and `Forcefield.create_prepped_molecule` 
    2. Forcefield parameters are stored in ParmEd objects instead of text files
    3. atom- and bond-specific terms can be retrieved through the `Atom.ff` and `Bond.ff` attributes
    4. Replaced `mdt.parameterize` with `mdt.create_ff_parameters`, which creates a Forcefield
       object with parameters for a passed molecule.
  - Molecular **alignments** using principal moments of inertia and bond-bond alignment
  - Better **constraint handling**, making it easier to add, remove, and clear geometric constraints
  - **Drug discovery energy models** - mmff94, mmff94s, and Ghemical - available through the `OpenBabelPotential` model ([\#111](https://github.com/autodesk/molecular-design-toolkit/issues/111))

##### INFRASTRUCTURE CHANGES
  - Simultaneous Python 2/3 support ([\#150](https://github.com/autodesk/molecular-design-toolkit/issues/150))
  - `moldesign` no longer requires `nbmolviz`, making for a much lighter-weight installation
    ([\#140](https://github.com/autodesk/molecular-design-toolkit/issues/140))
  - Users will be prompted to install the correct version of `nbmolviz` if necessary; that
    package now automatically installs itself in more cases
  - MDT can be configured to run external packages locally, if you'd prefer not to use docker
  - More robust molecular data structures make it harder (albeit still not THAT hard) to create
    an inconsistent topological state
  - More automation, runtime checks and notebook UI options to make sure sure that
    everything is installed correctly
  - "CCC" demo server removed. To automatically download and run dependencies like OpenBabel and
    NWChem, docker must be installed locally
  - Centralized handling of external software interactions via the `moldesign.compute.packages`
    module
    
##### BUGS
Test coverage has gone from <40% in the last release to 88% in the current one; 
the test and deploy pipeline is now fully automated; and tests have been added for a variety
of corner cases. All this testing exposed a veritable cornucopia of bugs, a panoply of
off-by-one-errors, typos, race conditions and more. These have all been fixed, leaving the
code on a much more stable footing moving forward.



### 0.7.3 - October 17, 2016

##### NEW MODELING FEATURES
 - [\#33](https://github.com/autodesk/molecular-design-toolkit/issues/33) - Add DFT w/ gradients; MP2, CASSCF, CASCI w/out gradients
 - Constrained minimizations w/ SHAKE and scipy's SLQSP
 - Transition dipoles and oscillator strengths
 - GAFF parameterizer for small molecules -- `params = mdt.parameterize(mol)`
 - AM1-BCC and Gasteiger partial charge calculators: `mdt.calc_am1_bcc_charges` and
    `mdt.calc_gasteiger_charges`
 - Add PDB database and biomolecular assembly support for mmCIF files
 - [\#72](https://github.com/autodesk/molecular-design-toolkit/issues/72) - Add `moldesign.guess_formal_charges` and `moldesign.add_missing_data`
 - Excited and multi-state property calculations with CAS methods
 - Rename `build_bdna` to `build_dna_helix` and give access to all NAB helix types


##### OTHER ENHANCEMENTS
 - [\#78](https://github.com/autodesk/molecular-design-toolkit/issues/78) - `moldesign` now imports much more quickly
 - Add `GAFF` energy model to automate small molecule parameterization
 - Change Example 2 to show an absorption spectrum calculation
 - Add Example 4 on protein MD with a small ligand
 - Add Example 5: on constrained minimization and enthalpic barriers
 - Add Tutorial 3: QM data analysis
 - Show changelog and version check in the `mdt.about()` (aka `mdt.configure`) widget
 - Change moldesign.tools and moldesign.helpers modules into more rationally organized subpackages
 - `mdt.set_dihedral` can be called with two atoms in the same way as `mdt.dihedral`
 - Explicit parameter created to store wavefunction guesses
 - Better access to density matrix in wavefunction objects
 - Improved parsing support for PDB and mmCIF files

##### BUGFIXES
 - [\#61](https://github.com/autodesk/molecular-design-toolkit/issues/61) - fixed a KeyError when parsing PDBs with metal centers or ions
 - [\#74](https://github.com/autodesk/molecular-design-toolkit/issues/74) - Add function to create PDB files with correct TER records (used for TLeap input)
 - Better handling of chains with non-standard residues
 - `mdt.add_hydrogens` no longer creates structures with separated residues
 - Fix sign of dihedral gradient
 - Charge quantities now mostly have the correct units


### 0.7.2 - July 26, 2016

##### NEW MODELING FEATURES
 - Add trajectory geometry analysis functions (`traj.dihedral`, `traj.angle`, etc.)
 - Can now calculate angles, dists, and dihedrals by name within a residue
    (`residue.distance('CA','CG')`)
 - Calculate dihedral angles using only two atoms defining the central bond (MDT will infer
    infer the other two in a consistent way)

##### CHANGES
 - Completed tutorials

##### BUGFIXES
 - [\#28](https://github.com/autodesk/molecular-design-toolkit/issues/28): Fixed a rounding error and logging problems with OpenMM trajectory snapshots
 - [\#21](https://github.com/autodesk/molecular-design-toolkit/issues/21): Better bond orders to structures in Amber files, which don't store them
 - [\#20](https://github.com/autodesk/molecular-design-toolkit/issues/20): Store OpenMM force vector correctly


### 0.7.1 - July 20, 2016

##### BUGFIXES
  - [\#4](https://github.com/autodesk/molecular-design-toolkit/issues/4): Use public demo CCC server by default
  - [\#3](https://github.com/autodesk/molecular-design-toolkit/issues/3): Fix `python -m moldesign intro`

#### 0.7.0 - July 15, 2016
 - Initial public release
