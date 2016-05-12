##Buckyball

### Who
Buckyball provides a python toolkit for computational chemists to build advanced workflows. Jupyter notebook allows researchers to easily share their workflows, allowing them to expose their research as sharable, reproducible software artifacts.

### What

 Buckyball is designed to solve pain points in 3D molecular modeling, allowing computational chemists with a working knowledge of python to build powerful workflows.
 
  * **Getting 3D structures:** file i/o for most common formats; PDB database queries; SMILES converter; IUPAC name converter DNA strand builder
  * **Manipulating 3D structures**: symmetrization methods; interactive internal coordinates modification; interactive selections
  * **Visualization**: 2D and 3D static viewers; visual force field setup; visual selection widgets 
  * **Modeling**: QM (HF, DFT, semi-empirical, configuration interaction); MM (Amber-type force fields); energy minimizations
  * **Simulation**: Generalized energy minimizations; NVE, NVT, and NPT dynamics
  * **Analysis**: Intuitive, information-rich data structures for atoms, residues, molecules, force fields, and electronic wavefunctions; plotable access to geometric and time-series data


### Where

 * Chemoinformatics: OpenBabel, RDKit, CDK
 * Modeling: *Buckyball*, Atomic Simulation Environment
 * Protein modeling: Rosetta
 * Bioinformatics: BioPython
 
 
 
### Abstract

Buckyball is a free, open source Python toolkit for molecular modeling, visualization, and cloud computing. Buckyball enables end-to-end computational chemistry workflows by combining a suite of intuitive python APIs for chemical modeling with a suite of web-based tools for visualization and distributed computing. It is specifically designed to take advantage of the popular IPython/Jupyter notebook environment: for example, in a single notebook, a user can download and view a PDB protein structure; perform a quantum chemical minimization on its small molecule ligand; assign MM force field parameters; launch and track a remote molecular dynamics simulation; and visualize the resulting trajectory. Users access this functionality through Buckyball's abstracted, object-oriented API, which relies on a core of well-established computational chemistry packages (including OpenMM, NWChem, OpenBabel, AmberTools, and 3DMol.js). To obviate the need for users to compile or deploy this additional software, the toolkit uses several modern, web-based technologies. First, software containerization is used to automatically deploy software to arbitrary computing resources. Second, visualizations are rendered in modern web browsers' 2D and 3D graphics engines, enabling highly portable molecular visualizations. Third, the Cloud Batch Scheduler (also introduced here) provides a simple, abstract interface to distributed computing resources. Finally, integration with the Jupyter notebook leverages its pre-existing strengths for computational exploration and sharable workflows. We believe that this unique combination of modern technologies provides a powerful, streamlined environment for structural chemistry research.