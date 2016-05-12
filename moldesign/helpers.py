import collections

import moldesign as mdt


def get_all_atoms(*objects):
    """ Given Atoms, AtomContainers, lists of Atoms, and lists of AtomContainers,
    return a flat list of all atoms contained therein.

    A given atom is only returned once, even if it's found more than once.

    Args:
        *objects (moldesign.Atom OR moldesign.AtomContainer OR List[moldesign.Atom] OR
            List[moldesign.AtomContainer]): objects to take atoms from

    """
    atoms = collections.OrderedDict()

    for obj in objects:
        if isinstance(obj, mdt.Atom):
            atoms[obj] = None
        elif hasattr(obj, 'atoms'):
            atoms.update((x,None) for x in obj.atoms)
        else:
            for item in obj:
                if isinstance(item, mdt.Atom):
                    atoms[item] = None
                elif hasattr(item, 'atoms'):
                    atoms.update((x, None) for x in item.atoms)

    return mdt.AtomList(atoms.iterkeys())





# def get_residues(obj, **queries):
#     """
#
#     Args:
#         obj ():
#         **queries ():
#
#     Returns:
#
#     """
#     for residue in obj.residues:
#         pass
#
