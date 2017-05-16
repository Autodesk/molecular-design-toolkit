## Molecule subpackage
The `moldesign.molecule` subpackage contains class definitions for molecular data,
including atomic, molecular, biomolecular, and trajectory data structures. 

### Naming and indexing conventions

#### Biomolecules
We run into a lot of cases where there's more than one way to name or index something. For instance, do we store a molecule's residues as a normal python list, indexed at 0? Or do we use the sequence numbers from the PDB file?

In the end, buckyball makes both available. Biologically relevant objects - atoms, chains, and residues, will each have a:
 * `name` - descriptive name, _unique_ in its enclosing structure
 * `index` - this object's list index in its molecule
 * `pdbname` - the object's name in the PDB entry
 * `pdbindex` - atom serial # or residue sequence number
 

| object/ attribute | `Atom`  | `Residue`  | Chain  | Molecule |
|---|---|---|---|---|
| `x.name` | PDB Atom name (e.g., `OE3`)  | PDB Name + PDB index (e.g., `ALA302`) | PDB Index (aka "chain identifier" - `A`, `B`, `C`, etc.)  | Filename |
| `x.index`| Index in `mol.atoms,` i.e., `atom==mol.atoms[atom.index]`  | Index in `mol.residues`, i.e. `res==mol.residues[res.index]`  | Index in `mol.chains`, i.e., `chain==mol.chain[chain.index]` | n/a|
| `x.pdbindex`| aka "Atom serial number" from pdb file (usually 1-based)| aka "Residue sequence number" from pdb file (residue's position in the primary structure, 1-based)  | "Chain identifier" from PDB file  | n/a|
| `x.pdbname` | Same as `atom.name` | PDB residue name (e.g., `ALA`)  | Same as `chain.name`, `chain.pdbindex` - `A`, `B`, `C`, etc. | 4-letter PDB code |
| `str(x)`| `Atom {name} (index {index}) in  residue {res.name}, chain {chain.name}, mol {mol.name}` | `Residue {resname} (index {index}) in chain {chain.name}, mol {mol.name}` | `Chain {chain.name} in mol {mol.name}`  | `Molecule {name} (N chains, M residues, L atoms)` |


Here's how the python attributes correspond to an entry in a PDB file:
```  
  chain.name == chain.pdbname     atom.x   atom.y  atom.z 
                       |             |       |       |   
ATOM     46   CB   ARG A  43      12.924  87.757  96.420  0.50 37.26           O
          |   |    |      |                                                    |
atom.pdbindex |    |    residue.pdbindex                                   atom.elem
              |    | 
              |   residue.resname == residue.pdbname
   atom.name == atom.pdbname
```

#### Small molecules
Small molecules can come from a variety of sources with a variety of different metadata available. If a given molecule is provided with PDB-type metadata, we'll name and index it according to the biomolecule conventions above. Otherwise, a 'placeholder' residue and chain will be created to hold the atoms. If the atom names are just the names of the elements (e.g., all carbon atoms are named C), atom.name will be automatically assigned as `"%s%d" % (atom.elem, atom.index)`.


### Data model specs

The MDT API is designed to let you manipulate molecules from a variety of levels, with internal components that keep the data consistent. For instance, we do the bookeeping to make sure that these statements are always true:
  - `residue.atoms[3].residue == residue`
  - `atom.position == mol.positions[atom.index]`
 
even if you manipulate the structure.
 
This means that any time the  user modifies data, the entire molecule needs to be self-consistently updated. `Molecule` objects are therefore tightly coupled with their `Atom`, `Residue`, and `Chain` components to achieve this.

Here are some general specifications and test for how the molecular data model reacts to changes:


##### List indices

Regardless of structure manipulation, it will always be true that

```python
atom is molecule.atoms[atom.index]
residue is molecule.residues[residue.index]
chain is molecule.chains[chain.index]
```

##### Adding a new object to an existing molecule
Either by settings its properties:
```
atom.molecule = mol
residue.molecule = mol [...]
```

Or equivalently:
```
mol.atoms.append(atom)
mol.residues.append(residue)
mol.chains.append(chain)
```

All of the object's children and parents (i.e., the associated chains, residues, atoms) will be assigned to this molecule if they're not already (`ValueError` is raised if any one of these objects is already different `Molecule`).

If `atom.residue is None`, the atom will made a child of the molecule's default residue (stored at `molecule._defresidue`). Similarly, unassigned residues become children of `molecule._defchain`.

 
##### Moving an atom between residues
An atom can be moved between residues by simple reassignment:
```python
atom = mol.residues[0].atoms[0]
atom.residue = mol.residues[3]
assert atom in mol.residues[3]
assert atom not in mol.residues[0]
```

##### Adding a new residue
Adding an atom to a new residue will automatically add that residue to the same molecule in a new chain:
```python
initmol = mol.copy()

newres = mdt.Residue(name='UNL')
atom.residue = newres

assert newres.molecule is atom.molecule
assert newres.chain.atoms == [atom]
assert mol.num_chains == initmol.num_chains + 1
assert mol.num_residues == initmol.num_residues + 1
assert mol.num_atoms == initmol.num_atoms
```

##### Changing an atom's name
Changing an atom's name will automatically update its residue:
```python
atom = mol.atoms[0]
res = atom.residue
oldname = atom.name

assert oldname in res
atom.name ='fakename'
assert atom.name in res
assert oldname not in res
assert res['fakename'] == atom
```


##### Removing an atom from a molecule
You can do this in two equivalent ways:
```python
atom.molecule = None
## OR ##
mol.atoms.remove(atom)
```
In both cases, the atom will automatically be removed from all primary and secondary structure, and the bonds it's part of will disappear as well.

However, empty residues and chains will remain in the molecule's primary structure even if they contain no atoms.


```
atom = mol.atoms[3]
oldres = atom.residue
oldchain = atom.chain
mol.atoms.remove(atom)

assert atom.residue is atom.chain is atom.molecule is None
assert atom not in oldres
assert atom not in oldchain
assert atom not in molecule
```

##### Removing residues and chains from a molecule
These objects are removed similarly to atoms:
```python
chain.molecule = None
## OR ##
mol.residues.remove(residue)
```
In both cases, these structures and all other structures they contain will be removed from the molecule as well.
```
residue = mol.residue[3]
mol.residues.remove(residue)

assert residue.chain is None

for atom in residue.atoms:
    assert atom not in mol
    assert atom in residue
    assert atom.chain is None
```


#### Behaviors
Assigning a atom to a molecule (`molecule.atoms.append(atom)` or `atom.molecule=mol`):
  - IF the atom is already part of this molecule
     - raise Exception (can't have atom in list twice)
  - IF the atom is part of a different molecule
     - raise Exception (can't be part of two different molecules)
  - IF the atom does NOT have a residue:
     - automatically assign it to the molecule's default residue
  - IF the atom is part of a residue and/or chain
     - raise Excpetion (need to assign the entire residue/chain)
  - IF the atom is bonded to atoms outside of this molecule
     - raise Exception (need to assign them all in one go)
     
Assigning a list of atoms to a molecule:
  - The order of the atoms will be preserved, even if it conflicts with the order of the residues
     
Removing an atom from a molecule:
  - It is removed from its parent residue
  - All bonds it is involved in will disappear from the molecule
  - However, the bonds are RETAINED by the unassigned atom
  - they bonds will re-appear if the atom is added back into the
    molecule, even in a different place
     
Creating a new bond:
  - If the atom is unassigned
     - it can be bonded to anything
     - its bonds will appear in a molecule when added to that moleclue.
  - If the atom belongs to a molecule
     - it can only be bonded to atoms in the same molecule
     - the bonds can cross residue/chain boundaries without restriciton
     
     
Assigning an atom a residue (`res.atoms.append(atom)` or `atom.residue=res`):
  - the atom automatically gets the same molecule as res (including `None`) 
  - the atom is removed from any previous molecule
  - the atom is added to residue's molecule (if present) at the end of its current residue

Adding a residue to a molecule (`residue.molecule=mol`) is same as with atoms, except:
 - IF the residue has a chain:
    - raise Exception (need to assign its entire chain)
 - ALL of the residues atoms will be assigned as well