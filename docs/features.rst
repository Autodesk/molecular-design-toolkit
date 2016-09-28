Molecular Design Toolkit Features
=================================

.. contents::
  :depth: 5



Molecular Mechanics
-------------------

 - :class:`models.Forcefield` Generic Forcefield energy model (automatically chosen implementation)

 - :class:`models.GAFF` Amber force field for small organic molecules

 - :class:`models.OpenMMEnergyModel` OpenMM energy kernel for Amber/CHARMM-like forcefields

 - :class:`models.SanderEnergyModel` Sander energy kernel (from AmberTools) for Amber/CHARMM-like forcefields

 - :class:`models.Spring` A single harmonic bond 2 atoms (for testing and toy models)

 - :class:`models.HarmonicOscillator` A 1-D Harmonic oscillator centered at x=0 (for testing and toy models)


Quantum Chemistry
-----------------

 - :class:`models.RHF` Restricted Hartree-Fock (automatically chosen implementation)

 - :class:`models.DFT` Density functional theories (automatically chosen implementation)

 - :class:`models.CASSCF` Complete-active space MCSCF (automatically chosen implementation)

 - :class:`models.Semiempirical` Generic semiempirical theories (automatically chosen implementation)

 - :class:`models.PySCF` PySCF *ab initio* QM energy kernel, includes implementations for a large number of quantum energy models, including RHF, DFT, CASCI, CASSCF, MP2, and coupled cluster

 - :class:`models.SQM` SQM semi-empirical energy kernel (from AmberTools), includes implementations for a large number of semiempirical models, including MNDO, AM1, PM3, PM6, DFTB


Molecular Dynamics
^^^^^^^^^^^^^^^^^^

 - :class:`integrators.VelocityVerlet` Generic Velocity Verlet dynamics (automatically chosen implementation)

 - :class:`integrators.Langevin` Generic Langevin type (automatically chosen implementation)

 - :class:`integrators.OpenMMVerlet` Velocity Verlet dynamics for use with OpenMM energy models

 - :class:`integrators.OpenMMLangevin` Velocity Langevin dynamics for use with OpenMM energy models

 - :class:`integrators.SurfaceHopping` Multi-state surface hopping dynamics using fewest switched. Implementation: internal.


Interactive visualization
-------------------------

 - :class:`viewer.Configurator` Automatically generates user interfaces for configuring simulations

 - :class:`viewer.GeometryViewer` 3D molecular viewer

 - :class:`viewer.ChemicalGraphViewer` 2D molecular viewer

 - :class:`widgets.OrbitalViewer` 3D molecular orbital viewer

 - :class:`widgets.BondSelector` widget for building lists of atoms and/or bonds

 - :class:`widgets.ResidueSelector` widget for building lists of atoms and/or residues

 - :class:`widgets.GeometryBuilder` widget for manipulating internal coordinates

 - :class:`widgets.Symmetrizer` widget for displaying and manipulating molecular symmetry groups

 - :class:`widgets.ParameterizationDisplay` 3D display of issues when assigning forcefield parameters


Data analysis
-------------

Simulation results are stored in numpy arrays with an explicit unit system based on ``pint`` for easy analysis and comparison. A few are shown here as examples:

**Static properties**

 - :meth:`Molecule.potential_energy`, :meth:`Molecule.forces`, :meth:`Molecule.dipole_moment`, :meth:`Molecule.wfn` Molecular properties calcualted by energy models: the potential energy, force array, dipole moment vector, and electronic wavefunction, respectively.

 - :class:`orbitals.ElectronWfn` A data structure storic electronic wavefunction information, as calculated by a quantum energy model.

 - :meth:`ElectronWfn.aobasis.fock` :meth:`ElectronWfn.aobasis.overlaps` The Fock and overlap matrices in the AO basis

 - :meth:`ElectronWfn.canonical.fock` :meth:`ElectronWfn.canonical.overlaps` The Fock and overlap matrices in the canonical orbital basis

 - :meth:`ElectronWfn.canonical.coeffs` The canonical orbital coefficients in the AO basis

**Trajectory properties**

 - :class:`Trajectory` A data structure storing a series of molecular structures with associated properties

 - :meth:`Trajectory.rmsd` Calculate a timeseries of RMSD values over the course of a trajectory

 - :meth:`Trajectory.distance`, :meth:`Trajectory.angle`, :meth:`Trajectory.dihedral` Return a timeseries of distances, angles, or dihedral angles over the course of a trajectory

 - :meth:`Trajectory.time` :meth:`Trajectory.potential_energy` :meth:`Trajectory.kinetic_temperature` ``...`` - Return timeseries of times, energies, temperatures, etc. over the course of a trajectory


Interfaces
----------

**Files and databases**

 - :meth:`read`, :meth:`write` read/write molecular file formats. Supports PDB, mmCIF, SDF, XYZ, MOL2, and pickled objects. Implementations: OpenBabel, BioPython, or internal.

 - :meth:`from_smiles` Convert a SMILES string into an MDT molecule with a 3D structure. Implementation: OpenBabel.

 - :meth:`from_name` Convert an IUPAC chemical name into an MDT molecule with a 3D structure. Implementation: Opsin.

 - :meth:`from_pdb` Download and create a molecule object from a PDB code. Implementation: BioPython.


**Python objects**

MDT molecules can also be converted into objects for a variety of other Python chemistry libraries:

 - :meth:`interfaces.mol_to_pybel`, :meth:`interfaces.pybel_to_mol` Convert an MDT molecule to/from a `pybel` (i.e. OpenBabel) molecule object.

 - :meth:`interfaces.mol_to_pyscf`, :meth:`interfaces.pyscf_to_mol` Convert an MDT molecule to/from a PySCF molecule object.

 - :meth:`interfaces.topology_to_mol`, :meth:`interfaces.mol_to_topology`  Convert an OpenMM topology object to/from an MDT molecule


Tools
-----

**Topology manipulation**

 - :meth:`add_hydrogen` Saturate a molecule's valence with hydrogens. Implementation: OpenBabel.

 - :meth:`guess_bond_orders` Assign bond orders based on geometry and/or topology. Implementation: OpenBabel.

 - :meth:`mutate_residues` Mutate DNA bases and amino acid residues. Implementation: PDBFixer.

 - :meth:`add_water_box` Add water box with optional ions. Implementation: PDBFixer.


**Forcefields**

 - :meth:`assign_forcefield` Returns a new molecule with forcefield assignments and any missing atoms. Implementation: AmberTools/tLeap.

 - :meth:`parameterize` Assign forcefield parameters to a molecule. Implementation: Ambertools/antechamber.

 - :meth:`calc_am1_bcc_charges` :meth:`calc_gasteiger_charges` :meth:`calc_esp_charges` Calculate partial charges for use with a forcefield. Implementation: Ambertools/antechamber and Ambertools/SQM
