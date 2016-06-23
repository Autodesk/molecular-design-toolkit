
The Molecule Object
-------------------
You'll spend most of your time in MDT working with :class:`Molecules <moldesign.Molecule>`. They
contain all the information necessary to describe a molecular system.

   >>> import moldesign as mdt
   >>> benzene = mdt.from_name('benzene')
   >>> protease = mdt.from_pdb('3AID')
   >>> type(benzene)
   mdt.Molecule
   >>> type(protease)
   mdt.Molecule

Below, we describe some of the most relevant attributes of a :class:`Molecule <moldesign.Molecule>` object.

Atoms
"""""
Each :class:`Molecule <moldesign.Atom>` contains a list of :class:`Atoms <moldesign.Atom>`.

   >>> len(benzene.atoms)
   12
   >>> atom = benzene.atoms[3]
   >>> atom.index
   3

Coordinates
"""""""""""
Coordinates for the entire :class:`Molecule <moldesign.Atom>`
are stored as flat vectors of length `3*mol.num_atoms`.

    >>> len(benzene.positions)
    36
    >>> len(benzene.momenta)
    36


Bonds
"""""
Bonds between atoms are stored as a graph, accessible via
:meth:`mol.bond_graph <moldesign.Molecule.bond_graph>`.

Primary structure
"""""""""""""""""
Biomolecules also contain primary structure information such as :class:`Chains <moldesign.Chain>`
and :class:`Residues <moldesign.Residue>`. :class:`Chains <moldesign.Chain>` can be accessed by name OR by index:

   >>> chain1 = protease.chains['A']
   >>> chain2 = protease.chains[0]
   >>> chain1 is chain2
   True

:class:`Residues <moldesign.Residue>` can similarly be accessed through a flat list or by name:

   >>> res0 = protease.residues[0]
   >>> resA = protease.chains['A'].residues['PRO1']
   >>> res0 is resA
   True

Molecular properties
""""""""""""""""""""
:class:`Molecules <moldesign.Molecule>` store :class`MolecularProperties` - groups of properties
that have been calculated by an :class:`EnergyModel` at the molecule's current position.

    >>> benzene.set_potential_model(mdt.models.RHF(basis='3-21g'))
    >>> benzene.calculate()
    >>> benzene.properties.potential_energy
    # [quantity with energy units]
    >>> len(benzene.properties.forces)
    36

Properties almost always include ``potential_energy``; other common properties include ``forces``, ``electronic_state``, and ``dipole``.

Note:
   These properties are only accessible if they correspond to the molecule's current position -
   they won't be returned if ``benzene.position`` changes.

Electronic structure
""""""""""""""""""""
Quantum chemical :class:`EnergyModels <moldesign.models.EnergyModelBase>` will also create an object representing the electronic wavefunction, accessible at :class:`mol.electronic_state (see its documentation for more details) <moldesign.orbitals.ElectronicWfn>`.

    >>> wfn = benzene.electronic_state
    >>> wfn.aobasis
    >>> wfn.molecular_orbitals['canonical']


