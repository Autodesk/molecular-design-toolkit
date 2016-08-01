# Programming guidelines

### Principles
 1. Users are **scientists, not programmers**. They'll need to know enough python to accomplish their science, but Molecular Design Toolkit (MDT) users may never know (or want to know) what a metaclass is or even how a web page works. Most workflows should be accomplishable without deep python or programming skills.
 1. Build **power tools for competent computational chemists** - users need customizable, introspectable, and composable tools to build complex simulation workflows. We can support users with good documentation, great examples, and sensible defaults (e.g., well-chosen convergence parameters, proper timesteps).
 1. Prefer a **clean, stable API** above all else. MDT's python API *is* its user interface. A good implementation is one that produces a good interface.
 1. When in doubt, strive for **usability and user-friendliness**

### Make information accessible and intuitive
MDT is designed, as much as possible, to allow users to interactively explore molecular systems. Here are some of the ways we try to make the API an intuitive experience:
 1. **Write for autocomplete**: users should be encouraged to type a couple letters then hit `tab`, even if they don't know exactly what they're looking for 
 1. **Keep namespaces flat**: don't make users hunt for a piece of data - make it accessible from the highest-level logical location (`mol.aobasis`, not  `mol.electronic_wavefunction.orbitals.aobasis`)
 1. **Add specific, use-case-based sugar** - e.g., users will frequently create lists of different types of residues. It's easier to type `waters = mol.get_residues(type='water')` than it is to type `waters = [res for res in mol.residues if res.type == 'water']`. If this is a *really common* use case, consider adding `mol.get_water_residues()`.
 1. **Use readable names**: long, descriptive names (`mol.calculate_potential_energy()`) are preferable to short, cryptic ones (`mol.calc_pe()`). Using jupyter means that everyone has access to autocomplete.
 1. **Give users lists, not iterators**: iterators are the enemies of easily-accessible data.  No one wants to type `list(mol.atoms())[3]`, or even worse, `[atom for atom in mol.atoms() if atom.idx==3][0]`. Expose an explicit iterator method if you really need to.
 1. **Expose lightweight data as attributes and with `@property`** - Users should be able to access any available data - the number of atoms, precalculated forces, or a molecule's kinetic energy, for instance - in a single line of code *without an explicit function call*. For many molecular quantities (such as the number of degrees of freedom in a molecule), we use descriptors and `@property`s to make them immediately accessible to the user. The principle data structures - `Molecule`s, `Atom`s and `BioUnit`s - all make liberal use of descriptors and `@property`s.
     * NO: `ke = mol.get_kinetic_energy()`
     * YES: `ke = mol.kinetic_energy`
     * NO: `chain = list(mol.getChains())[0]; residue = chain.getResidues(seqidx=35)[0]`
     * YES: `res = mol.chains['A'].residues['ALA35']`
 1. **Expose heavy-duty computation as methods** - A potential energy calculation could take anywhere from milliseconds to days. This should NOT be triggered by a user accessing `mol.potential_energy`; instead, they should call `mol.calculate_potential_energy()`.

### Keeping the API clean, flat, and Jupyter-friendly
1. The **top level `moldesign` namespace** should contain everything the user will need:<br> 
   - YES: `moldesign.[name]`, e.g. `moldesign.from_smiles`
   - NO: `moldesign.[name1].[name2].[name3]`, e.g. `moldesign.interfaces.openbabel.from_smiles`
1. Use **delegation** to flatten complex objects: override `__getattr__` to allow an object to call its attributes' methods.<br>
   - YES: `trajectory.set_style('vdw')`
   - NO: `trajectory.viewer.set_style('vdw')`
1. Make sure **autocomplete and pop-up docstrings** are useful:
   - If you override `__getattr__`, also override `__dir__` to enable Jupyter's autocomplete functions.
   - Use **explicit function signatures** whenever possible (e.g., avoid `my_func(*args, **kwargs)`) so that Jupyter can provide a useful call signature. (possibly use the `decorators` module to copy call signatures?)
1. **Synonyms** - users shouldn't have to look at the documentation to figure out if it's `mol.numatoms`, `mol.natoms`, or `mol.num_atoms`. It's totally fine, even preferable, to have these all be the same thing (use `@property` to make sure that they are all equivalent).
1. **Keep links to related objects** - if an atom belongs to a molecule, make it accessible as `atom.molecule`. This will create lots of circular composition (cases where `object1.child is child.object1`) - this is fine.
1. **Don't clutter namespaces**:
   - *objects:* use underscores to hide objects' internal variables:<br>*ex:* `trajectory._position_cache`, not `trajectory.position_cache`
   - *modules:* use the `__all__` attribute in all submodules to define what gets imported to the top-level namespace.
 
### Prefer composition to inheritance
Even advanced programmers have difficulty understanding complex inheritance trees, and python's syntax for passing arguments to superclasses only adds to the confusion. Whenever possible, build complex objects using composition instead of inheritance. Use attribute delegation (override `__getattr__` if necessary - see above) to keep the object's namespace flat.

For example, the `TrajectoryViewer` class incorporates both a `GeometryViewer` and several ipywidgets. These are all set up as attributes of the `TrajectoryViewer` class: 
```python
trajview = TrajectoryViewer( mol )
trajview.viewer  # GeometryViewer object
trajview.slider  # ipywidget float slider
trajview.frameinspector  # bb.ui.FrameInspector object
```

To keep the namespace flat, `TrajectoryViewer` objects delegate calls to their viewer: `TrajectoryViewer.vdw()` is a synonym for `TrajectectoryViewer.viewer.vdw()`.<br>Similarly, for a `Molecule` object `mol = Molecule()`,  `Molecule.atomic_number` is a synonym for `[atom.atomic_number for atom in Molecule.atoms]`.

### Some inheritance is OK
Use inheritance when it makes sense and is easy to understand. For example:

1. Use **abstract base classes** to define interfaces -- e.g., all energy models derive from the abstract `EnergyModelBase`, because they all offer `.prep`, `.calculate`, `.DEFAULT_PROPERTIES`, etc. Similarly, `Residue` and `Chain` both inherit from `BioUnit`.
   * But don't use the built-in `abc` module, which requires too much conceptual overhead for too little benefit
1. **Mixins** - defining a set of methods that should work with a variety of classes (e.g., the `AtomContainer` class is mixed into `Atom`, `AtomList`, `BioUnit`, and `Molecule`, giving all of these classes access to `self.distance`, `self.copy_topology`. We're also using "trivial mix-ins" (designed just to mix together one big class) to help organize the `Atom` and `Molecule` types.
1. **When it makes sense** - complex inheritance is preferable to code duplication. Sometimes it really is the best answer.
 
### Code style
1. Functions and variables should be `lowercase_with_underscores`. Class names and constructors should be `CapitalizedCamelCase`.
1. The user-exposed API should be clean, [PEP 8](https://www.python.org/dev/peps/pep-0008/) code.
1. Internally, readability and functionality are more important than consistency - that said,  [Google's style guide](https://google.github.io/styleguide/pyguide.html), along with [PEP 8](https://www.python.org/dev/peps/pep-0008/), is strongly encouraged.


# Contributing

### Who should contribute?
Anyone with a molecular modeling workflow that they want to enable or share. Experience and research-level knowledge of the field is an important asset! In contrast, limited programming experience *is definitely not* barrier to contributing - we can help you with that! Please ask for help getting started in our forums [link].

### Scope: What goes into MDT?
Established techniques and general simulation tools that will be useful for **3-dimensional biomolecular modeling**. MDT aims to enable scientists to easily build new simulation techniques and workflows, but new, immature techniques, or those with limited applicability outside of a particular system should be implemented as separate projects that *use* MDT, not *part of* MDT.
###### Could (and should!) be implemented in MDT:
 * Physical simulation and modelilng: Lambda dynamics; homology modelling; surface hopping; RPMD; metadynamics; markov state models; a library of common 3D structures (such as amino acids, carbon nanotubes, small molecules, etc.)
 * Visualization and UI: transitions between different views; interactive structure building and editing; ray-traced rendering; movie exports
 
###### Should implemented as a separate project:
 * Computational techniques: fluid dynamics solver (not useful at the atomic level), biological network models (no clear connection to 3D structures); machine-learning based quantum chemistry (immature, untested)
 * Visualization and UI: visualizations for specific systems (not generally applicable); 


