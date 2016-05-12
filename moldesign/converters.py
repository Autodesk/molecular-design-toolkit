import os
import cStringIO as StringIO
import cPickle as pkl
import gzip, bz2
import string
import functools

import moldesign as mdt
import moldesign.interfaces.openbabel as obabel
from moldesign.interfaces.opsin_interface import name_to_smiles
from moldesign.interfaces.openmm import amber_to_mol as read_amber
import moldesign.interfaces.biopython_interface as biopy
from moldesign.interfaces.ambertools import build_bdna

__all__ = ('read from_smiles from_pdb from_name build_bdna write read_amber '
           'build_assembly').split()

# TODO: if cpickle fails, retry with regular pickle to get a better traceback
PICKLE_EXTENSIONS = set("p pkl pickle mdt".split())
COMPRESSION = {'gz': gzip.open,
               'gzip': gzip.open,
               'bz2': bz2.BZ2File,
               'bzip2': bz2.BZ2File}  # extensions must be lower case here
READERS = {'pdb': biopy.parse_pdb,
           'cif': biopy.parse_pdb}
WRITERS = {}
for ext in PICKLE_EXTENSIONS:
    READERS[ext] = pkl.load
    WRITERS[ext] = functools.partial(pkl.dump, protocol=2)

from_smiles = obabel.from_smiles


def read(f, format=None):
    """ Read in a molecule from a file, file-like object, or string.
    Will also depickle a pickled object.

    Note:
        Files with ``.bz2`` or ``.gz`` suffixes will be automatically decompressed.
        Currently does not support files with more than one record - only returns the first record

    Args:
        f (str or file-like): Either a path to a file, OR a string with the
            file's contents, OR a file-like object
        format (str): molecule format (pdb, xyz, sdf, etc.) or pickle format
            (recognizes p, pkl, or pickle); guessed from filename if not passed

    Returns:
        moldesign.Molecule or object: molecule parsed from the file (or python object, for
            pickle files)
            
    Raises:
        ValueError: if ``f`` isn't recognized as a string, file path, or file-like object
    """
    filename = suffix = filestring = fileobj = None

    # Figure what sort of object was passed
    if isinstance(f, str) and os.path.exists(f):
        filename = f
        suffix = filename.split('.')[-1]
        if suffix in COMPRESSION:  # deal with compressed file paths
            fileobj = COMPRESSION[suffix](filename, 'r')
            suffix = filename.split('.')[-2]
        if format is None:
            format = suffix
    elif hasattr(f, 'read'):
        fileobj = f
    elif isinstance(f, str):
        filestring = f

    # Open a file-like object for the file
    if fileobj is None:
        if filestring is not None:
            fileobj = StringIO.StringIO(filestring)
        elif filename is not None:
            fileobj = open(filename, 'r')
        else:
            raise ValueError('Parameter to moldesign.read (%s) not ' % str(f) +
                             'recognized as string, file path, or file-like object')

    if format in READERS:
        mol = READERS[format](fileobj)
    else:  # default to openbabel if there's not an explicit writer for this format
        mol = obabel.read_stream(fileobj, format)

    if filename is not None:
        mol.name = filename

    return mol


def write(obj, filename=None, format=None, mode='w'):
    """ Write a molecule to a file or string.
    Will also pickle arbitrary python objects.

    Note:
        Files with ``.bz2`` or ``.gz`` suffixes will be automatically compressed.

    Args:
        obj (moldesign.Molecule or object): the molecule to be written
            (or python object to be pickled)
        filename (str): filename (if not passed, then a string is returned)
        format (str): molecule format (pdb, xyz, sdf, etc.) or a pickle file extension
            ('pkl' and 'mdt' are both accepted)

    Returns:
        str: if filename is none, return the output file as a string (otherwise returns ``None``)
    """

    # create the file object to write to
    if filename:
        return_string = False
        suffix = filename.split('.')[-1]
        if suffix in COMPRESSION:
            fileobj = COMPRESSION[suffix](filename, mode=mode)
            suffix = filename.split('.')[-2]
        else:
            fileobj = open(filename, mode=mode)
        if format is None:
            format = suffix
    else:
        return_string = True
        fileobj = StringIO.StringIO()

    # write to fileobj
    if format in WRITERS:
        WRITERS[format](obj, fileobj)
    else:
        fileobj.write(obabel.write_string(obj, format))

    if return_string:
        return fileobj.getvalue()
    else:
        fileobj.close()


def from_pdb(pdbcode):
    """ Import the given molecular geometry from PDB.org
        
    See Also:
        http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction

    Args:
        pdbcode (str): 4-character PDB code (e.g. 3AID, 1BNA, etc.)

    Returns:
        moldesign.Molecule: molecule object
    """
    import urllib2
    assert len(pdbcode) == 4, "%s is not a valid pdb ID." % pdbcode
    request = urllib2.urlopen('http://www.rcsb.org/pdb/files/%s.pdb' % pdbcode)
    ":type: urllib2.req"
    filestring = request.read()
    mol = read(filestring, format='pdb')
    mol.name = pdbcode
    return mol


def from_name(name):
    """Attempt to convert an IUPAC or common name to a molecular geometry.

    Args:
        name (str): molecular name (generally IUPAC - some common names are also recognized)

    Returns:
        moldesign.Molecule: molecule object
    """
    # TODO: fallback to http://cactus.nci.nih.gov/chemical/structure
    smi = name_to_smiles(name)
    mol = from_smiles(smi, name)
    return mol


def build_assembly(mol, assembly_name):
    """ Create biological assembly using a bioassembly specification.

    This routine builds a biomolecular assembly using the specification from a PDB header (if
    present, this data can be found in the  "REMARK 350" lines in the PDB file). Assemblies are
    author-assigned structures created by copying, translating, and rotating a subset of the
    chains in the PDB file.

    See Also:
        http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/biological-assemblies
        
    Args:
        mol (moldesign.Molecule): Molecule with assembly data (assembly data will be created by the
            PDB parser at ``molecule.properties.bioassembly``)
        assembly_name (str): name of the biomolecular assembly to build.

    Returns:
        mol (moldesign.Molecule): molecule containing the complete assembly

    Raises:
        AttributeError: If the molecule does not contain any biomolecular assembly data
        KeyError: If the specified assembly is not present
    """
    if 'bioassemblies' not in mol.properties:
        raise AttributeError('This molecule does not contain any biomolecular assembly data')
    try:
        asm = mol.properties.bioassemblies[assembly_name]
    except KeyError:
        raise KeyError(('The specified assembly name ("%s") was not found. The following '
                        'assemblies are present: %s') %
                       (assembly_name,
                        ', '.join(mol.properties.bioassemblies.keys())))

    # Make sure each chain gets a unique name - up to all the letters in the alphabet, anyway
    used_chain_names = set()
    alpha = iter(string.ascii_uppercase)

    # Create the new molecule by copying, transforming, and renaming the original chains
    all_atoms = mdt.AtomList()
    for i, t in enumerate(asm.transforms):
        for chain_name in asm.chains:
            chain = mol.chains[chain_name].copy()
            chain.transform(t)

            while chain.name in used_chain_names:
                chain.name = alpha.next()
            used_chain_names.add(chain.name)
            chain.pdbname = chain.pdbindex = chain.name
            all_atoms.extend(chain.atoms)
    newmol = mdt.Molecule(all_atoms,
                          name="%s (bioassembly %s)" % (mol.name, assembly_name))
    return newmol