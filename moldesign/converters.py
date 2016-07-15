# Copyright 2016 Autodesk Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
import os
import cStringIO as StringIO
import cPickle as pkl # TODO: if cpickle fails, retry with regular pickle to get a better traceback
import gzip
import bz2
import string
import functools

import moldesign as mdt
import moldesign.interfaces.openbabel as obabel
import moldesign.molecules.atomcollections
from moldesign.interfaces.opsin_interface import name_to_smiles
import moldesign.interfaces.biopython_interface as biopy

# These routines are imported here solely to make them available in this namespace
from moldesign.interfaces.openbabel import mol_to_pybel
from moldesign.interfaces.openmm import amber_to_mol as read_amber
from moldesign.interfaces.ambertools import build_bdna


def exports(o, name=None):
    __all__.append(o.__name__)
    return o
__all__ = ['read_amber', 'mol_to_pybel', 'build_bdna', 'from_smiles']



# NOTE - all extensions MUST be lower case here
PICKLE_EXTENSIONS = set("p pkl pickle mdt".split())
COMPRESSION = {'gz': gzip.open,
               'gzip': gzip.open,
               'bz2': bz2.BZ2File,
               'bzip2': bz2.BZ2File,
               None: open}
READERS = {'pdb': biopy.parse_pdb,
           'cif': biopy.parse_pdb}
WRITERS = {}
for ext in PICKLE_EXTENSIONS:
    READERS[ext] = pkl.load
    WRITERS[ext] = functools.partial(pkl.dump, protocol=2)

from_smiles = obabel.from_smiles


@exports
def read(f, format=None):
    """Read in a molecule from a file, file-like object, or string.
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
    filename = None

    # Open a file-like object
    if isinstance(f, str) and os.path.exists(f):  # it's a path to a file
        filename = f
        format, compression = _get_format(filename, format)
        fileobj = COMPRESSION[compression](filename, mode='r')
    elif hasattr(f, 'open'):  # we can get a file-like object
        fileobj = f.open('r')
    elif hasattr(f, 'read'):  # it's already file-like
        fileobj = f
    elif isinstance(f, str):  # it's just a string
        fileobj = StringIO.StringIO(f)
    else:
        raise ValueError('Parameter to moldesign.read (%s) not ' % str(f) +
                         'recognized as string, file path, or file-like object')

    if format in READERS:
        mol = READERS[format](fileobj)
    else:  # default to openbabel if there's not an explicit reader for this format
        mol = obabel.read_stream(fileobj, format)

    if filename is not None:
        mol.name = filename

    return mol


@exports
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
    # TODO: handle writing and return file-like objects
    format, compression = _get_format(filename, format)

    # First, create an object to write to (either file handle or file-like buffer)
    if filename:
        return_string = False
        fileobj = COMPRESSION[compression](filename, mode=mode)
    else:
        return_string = True
        fileobj = StringIO.StringIO()

    # Now, write to the object
    if format in WRITERS:
        WRITERS[format](obj, fileobj)
    else:
        fileobj.write(obabel.write_string(obj, format))

    # Return a string if necessary
    if return_string:
        return fileobj.getvalue()
    else:
        fileobj.close()


@exports
def write_trajectory(traj, filename=None, format=None, overwrite=True):
    """ Write trajectory a file (if filename provided) or file-like buffer

    Args:
        traj (moldesign.molecules.Trajectory): trajectory to write
        filename (str): name of file (return a file-like object if not passed)
        format (str): file format (guessed from filename if None)
        overwrite (bool): overwrite filename if it exists

    Returns:
        StringIO: file-like object (only if filename not passed)
    """
    format, compression = _get_format(filename, format)

    # If user is requesting a pickle, just dump the whole thing now and return
    if format.lower() in PICKLE_EXTENSIONS:
        write(traj, filename=filename, format=format)


    # for traditional molecular file formats, write the frames one after another
    else:
        tempmol = traj._tempmol
        if filename and (not overwrite) and os.path.exists(filename):
            raise IOError('%s exists' % filename)
        if not filename:
            fileobj = StringIO.StringIO()
        else:
            fileobj = open(filename, 'w')

        for frame in traj.frames:
            traj.apply_frame(frame)
            fileobj.write(tempmol.write(format=format))

        if filename is None:
            fileobj.seek(0)
            return fileobj
        else:
            fileobj.close()



@exports
def mol_to_openmm_sim(mol):
    try:
        return mol.energy_model.get_openmm_simulation()
    except AttributeError:
        raise AttributeError("Can't create an OpenMM object - no OpenMM energy_model present")


@exports
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


@exports
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


@exports
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
    all_atoms = moldesign.molecules.atomcollections.AtomList()
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


def _get_format(filename, format):
    """ Determine the requested file format and optional compression library

    Args:
        filename (str or None): requested filename, if present
        format (str or None): requested format, if present

    Returns:
        (str, str): (file format, compression format or ``None`` for no compression)

    Examples:
        >>> _get_format('mymol.pdb', None)
        ('pdb', None)
        >>> _get_format('smallmol.xyz.bz2', None)
        ('xyz','bz2')
        >>> _get_format('mymol.t.gz', 'sdf')
        ('sdf','gz')
    """
    if filename is None and format is None:
        raise ValueError('No filename or file format specified')
    elif filename is None:
        return format, None

    fname, extn = os.path.splitext(filename)
    suffix = extn[1:].lower()
    compressor = None
    if suffix in COMPRESSION:
        compressor = suffix
        suffix = os.path.splitext(fname)[1][1:].lower()

    if format is None:
        format = suffix

    return format, compressor
