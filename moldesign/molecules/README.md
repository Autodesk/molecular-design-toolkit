## Molecule subpackage
The `moldesign.molecule` subpackage contains class definitions for molecular data,
including atomic, molecular, biomolecular, and trajectory data structures. 

### Naming and indexing conventions

#### Biomolecules
We run into a lot of cases where there's more than one way to name or index something. For instance, do we store a molecule's residues as a normal python list, indexed at 0? Or do we use the sequence numbers from the PDB file?

In the end, buckyball makes both available. Biologically relevant objects - atoms, chains, and residues, will each have a:
 * `name` - a short (1-5 character) string
 * `index` - the index in python lists
 * `pdbname` - the official name from the PDB
 * `pdbindex` - the serial # (for atoms), sequence number (for residues), or chain identifier
 

| object/ attribute | `Atom`  | `Residue`  | Chain  | Molecule |
|---|---|---|---|---|
| `x.name` | PDB Atom name (e.g., `OE3`)  | PDB Name + PDB index (e.g., `ALA302`) | PDB Index (aka "chain identifier" - `A`, `B`, `C`, etc.)  | Filename |
| `x.index`| Index in `mol.atoms,` i.e., `atom==mol.atoms[atom.index]`  | Index in `mol.residues`, i.e. `res==mol.residues[res.index]`  | Index in `mol.chains`, i.e., `chain==mol.chain[chain.index]` | n/a|
| `x.pdbindex`| aka "Atom serial number" from pdb file (usually 1-based)| aka "Residue sequence number" from pdb file (residue's position in the primary structure, 1-based)  | "Chain identifier" from PDB file  | n/a|
| `x.pdbname` | Same as `atom.name` | PDB residue name (e.g., `ALA`)  | Same as `chain.name`, `chain.pdbindex` - `A`, `B`, `C`, etc. | 4-letter PDB code |
| `str(x)`| `Atom {name} (index {index}) in  residue {res.name}, chain {chain.name}, mol {mol.name}` | `Residue {resname} (index {index}) in chain {chain.name}, mol {mol.name}` | `Chain {chain.name} in mol {mol.name}`  | `Molecule {name} (N chains, M residues, L atoms)` |


Here's how the python attributes correspond to an entry in a PDB file:
```
ATOM     46   CB   ARG A  43      12.924  87.757  96.420  0.50 37.26           O
          |   |    |  |   |                                                    |
atom.pdbindex |    |  |  residue.pdbindex                                 atom.elem
              |    | chain.name == chain.pdbindex == chain.pdbname
              |   residue.pdbname
   atom.name == atom.pdbname
```

#### Small molecules
Small molecules can come from a variety of sources with a variety of different metadata available. If a given molecule is provided with PDB-type metadata, we'll name and index it according to the biomolecule conventions above. Otherwise, a 'placeholder' residue and chain will be created to hold the atoms. If the atom names are just the names of the elements (e.g., all carbon atoms are named C), atom.name will be automatically assigned as `"%s%d" % (atom.elem, atom.index)`.