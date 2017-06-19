from __future__ import print_function, absolute_import, division
from future.builtins import *
from future import standard_library
standard_library.install_aliases()

# Copyright 2017 Autodesk Inc.
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
from past.builtins import basestring
from future.utils import PY2

import bz2
import pickle as pkl  # TODO: if cpickle fails, retry with regular pickle to get a better traceback
import io
import functools
import gzip
import os

import moldesign as mdt
from . import utils
from .interfaces import biopython_interface
from .interfaces import openbabel as openbabel_interface
from .interfaces.parmed_interface import write_pdb, write_mmcif
from .helpers import pdb
from .external import pathlib

# imported names
read_amber = mdt.interfaces.openmm.amber_to_mol
from_smiles = mdt.interfaces.openbabel.from_smiles
from_inchi = mdt.interfaces.openbabel.from_inchi


utils.exports_names('from_smiles', 'read_amber', 'from_inchi')


@utils.exports
def read(f, format=None):
    """Read in a molecule from a file, file-like object, or string.
    Will also depickle a pickled object.

    Note:
        Files with ``.bz2`` or ``.gz`` suffixes will be automatically decompressed.
        Currently does not support files with more than one record - only returns the first record

    Args:
        f (str or file-like or pathlib.Path): Either a path to a file, OR a string with the
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
    closehandle = False
    streamtype = io.StringIO
    readmode = 'r'
    if PY2:
        f = pathlib._backport_pathlib_fixup(f)

    try:  # Gets a file-like object, depending on exactly what was passed:

        # it's a path to a file
        if _ispath(f):
            path = pathlib.Path(f).expanduser()
            format, compression, modesuffix = _get_format(path.name, format)
            fileobj = COMPRESSION[compression](str(path), mode='r' + modesuffix)
            closehandle = True

        # it can create a file-like object
        elif hasattr(f, 'open'):
            if format in PICKLE_EXTENSIONS:
                readmode = 'rb'
            fileobj = f.open(readmode)
            closehandle = True

        # it's already file-like
        elif hasattr(f, 'read'):
            fileobj = f

        # It's a string with a file's content
        elif isinstance(f, basestring):
            if format is None:
                raise IOError(('No file named "%s"; ' % f[:50]) +
                              'please set the `format` argument if you want to parse the string')
            elif format in PICKLE_EXTENSIONS or isinstance(f, bytes):
                streamtype = io.BytesIO
            fileobj = streamtype(f)

        else:
            raise ValueError('Parameter to moldesign.read (%s) not ' % str(f) +
                             'recognized as string, file path, or file-like object')

        if format in READERS:
            mol = READERS[format](fileobj)
        else:  # default to openbabel if there's not an explicit reader for this format
            mol = openbabel_interface.read_stream(fileobj, format)
    finally:
        if closehandle:
            fileobj.close()

    if filename is not None and mol.name not in (None, 'untitled'):
        mol.name = filename

    mdt.helpers.atom_name_check(mol)
    return mol


def _ispath(f):
    if not f:
        return False
    elif isinstance(f, basestring) and _isfile(f):
        return True
    elif isinstance(f, pathlib.Path):
        return True
    elif PY2 and pathlib._backportpathlib and isinstance(f, pathlib._backportpathlib.Path):
        return True
    else:
        return False


def _isfile(f):
    try:
        return os.path.isfile(f)
    except (TypeError, ValueError):
        return False


@utils.exports
def write(obj, filename=None, format=None):
    """ Write a molecule to a file or string.
    Will also pickle arbitrary python objects.

    Notes:
        Files with ``.bz2`` or ``.gz`` suffixes will be automatically compressed.
        Users will usually call this by way of ``Molecule.write``.

    Args:
        obj (moldesign.Molecule or object): the molecule to be written
            (or python object to be pickled)
        filename (str or pathlib.Path or file-like): path or buffer to write to (if not passed),
            string is returned
        format (str): molecule format (pdb, xyz, sdf, etc.) or a pickle file extension
            ('pkl' and 'mdt' are both accepted)

    Returns:
        str: if filename is none, return the output file as a string (otherwise returns ``None``)
    """
    close_handle = False

    if hasattr(filename, 'write'):  # it's a stream
        fileobj = filename
        return_string = False

        if format is None and getattr(fileobj, 'name', None):
            format, compression, mode = _get_format(fileobj.name, format)

        if format is None:
            raise ValueError("Could not determine format to write - please the 'format' argument")

    else:
        if ( format is None
             and filename is not None
             and not isinstance(filename, pathlib.Path)
             and len(str(filename)) < 5
             and '.' not in filename):
            # lets users call mdt.write(obj, 'pdb') and get a string (without needing the "format" keyword
            path, format = None, filename
        elif filename is None:
            path = None
        elif isinstance(filename, str):
            path = pathlib.Path(filename).expanduser()
        else:
            path = filename
            filename = str(path)

        writemode = 'w'
        format, compression, mode = _get_format(filename, format)

        if mode == 'b':
            streamtype = io.BytesIO
        else:
            streamtype = io.StringIO

        # First, create an object to write to (either file handle or file-like buffer)
        if path:
            return_string = False
            fileobj = COMPRESSION[compression](str(path), mode=writemode + mode)
            close_handle = True
        else:
            return_string = True
            fileobj = streamtype()

    # Now, write to the object
    if format in WRITERS:
        WRITERS[format](obj, fileobj)
    else:
        fileobj.write(openbabel_interface.write_string(obj, format))

    # Return a string if necessary
    if return_string:
        return fileobj.getvalue()
    elif close_handle:
        fileobj.close()


@utils.exports
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
    format, compression, modesuffix = _get_format(filename, format)

    # If user is requesting a pickle, just dump the whole thing now and return
    if format.lower() in PICKLE_EXTENSIONS:
        write(traj, filename=filename, format=format)

    # for traditional molecular file formats, write the frames one after another
    else:
        if filename and (not overwrite) and _isfile(filename):
            raise IOError('%s exists' % filename)
        if not filename:
            fileobj = io.StringIO()
        else:
            fileobj = open(filename, 'w' + modesuffix)

        for frame in traj.frames:
            fileobj.write(frame.write(format=format))

        if filename is None:
            fileobj.seek(0)
            return fileobj
        else:
            fileobj.close()


def read_pdb(f, assign_ccd_bonds=True):
    """ Read a PDB file and return a molecule.

    This uses ParmEd's parser to get the molecular structure, with additional functionality
    to assign Chemical Component Dictionary bonds, detect missing residues, and find
    biomolecular assembly information.

    Note:
        Users won't typically use this routine; instead, they'll use ``moldesign.read``, which will
        delegate to this routine when appropriate.

    Args:
        f (filelike): filelike object giving access to the PDB file (must implement readline+seek)
        assign_ccd_bonds (bool): Use the PDB Chemical Component Dictionary (CCD) to create bond
            topology (note that bonds from CONECT records will always be created as well)

    Returns:
        moldesign.Molecule: the parsed molecule
    """
    assemblies = pdb.get_pdb_assemblies(f)
    f.seek(0)
    mol = mdt.interfaces.parmed_interface.read_pdb(f)
    mol.properties.bioassemblies = assemblies

    f.seek(0)
    mol.metadata.missing_residues = mdt.helpers.get_pdb_missing_residues(f)

    # Assign bonds from residue templates
    if assign_ccd_bonds:
        pdb.assign_biopolymer_bonds(mol)

    if assemblies:
        pdb.warn_assemblies(mol, assemblies)

    return mol


def read_mmcif(f):
    """ Read an mmCIF file and return a molecule.

    This uses OpenBabel's basic structure parser along with biopython's mmCIF bioassembly parser

    Note:
        Users won't typically use this routine; instead, they'll use ``moldesign.read``, which will
        delegate to this routine when appropriate.

    Args:
        f (filelike): file-like object that accesses the mmCIF file (must implement seek)

    Returns:
        moldesign.Molecule: the parsed molecular structure
    """
    mol = mdt.interfaces.parmed_interface.read_mmcif(f)
    f.seek(0)
    assemblies = biopython_interface.get_mmcif_assemblies(f)
    if assemblies:
        pdb.warn_assemblies(mol, assemblies)
    mol.properties.bioassemblies = assemblies
    return mol


def read_xyz(f):
    tempmol = openbabel_interface.read_stream(f, 'xyz')
    for atom in tempmol.atoms:
        atom.residue = None
        atom.chain = None
    return mdt.Molecule(tempmol.atoms)


def write_xyz(mol, fileobj):
    fileobj.write(u"   %d\n%s\n" % (mol.num_atoms, mol.name))
    for atom in mol.atoms:
        x, y, z = atom.position.value_in(mdt.units.angstrom)
        fileobj.write(u"%s   %24.14f   %24.14f   %24.14f\n" % (atom.element, x, y, z))


@utils.exports
def from_pdb(pdbcode, usecif=False):
    """ Import the given molecular geometry from PDB.org
        
    See Also:
        http://pdb101.rcsb.org/learn/guide-to-understanding-pdb-data/introduction

    Args:
        pdbcode (str): 4-character PDB code (e.g. 3AID, 1BNA, etc.)
        usecif (bool): If False (the default), use the PDB-formatted file (default).
           If True, use the mmCIF-format file from RCSB.org.

    Returns:
        moldesign.Molecule: molecule object
    """
    import requests
    assert len(pdbcode) == 4, "%s is not a valid PDB ID." % pdbcode

    fileext = 'cif' if usecif else 'pdb'

    url = 'https://www.rcsb.org/pdb/files/%s.%s' % (pdbcode, fileext)
    request = requests.get(url)

    if request.status_code == 404 and not usecif:  # if not found, try the cif-format version
        print('WARNING: %s.pdb not found in rcsb.org database. Trying %s.cif...' % (
            pdbcode, pdbcode), end=' ')
        retval = from_pdb(pdbcode, usecif=True)
        print('success.')
        return retval

    elif request.status_code != 200:
        raise ValueError('Failed to download %s.%s from rcsb.org: %s %s' % (
            pdbcode, fileext, request.status_code, request.reason))

    filestring = request.text

    mol = read(filestring, format=fileext)
    mol.name = pdbcode
    mol.metadata.sourceurl = url
    mol.metadata.structureurl = 'https://www.rcsb.org/pdb/explore.do?structureId=%s' % pdbcode
    mol.metadata.pdbid = pdbcode
    if usecif:
        mol.metadata.sourceformat = 'mmcif'
    else:
        mol.metadata.sourceformat = 'pdb'
    return mol


@utils.exports
def from_name(name):
    """Attempt to convert an IUPAC or common name to a molecular geometry.

    Args:
        name (str): molecular name (generally IUPAC - some common names are also recognized)

    Returns:
        moldesign.Molecule: molecule object
    """
    from moldesign.interfaces.opsin_interface import name_to_smiles

    # TODO: fallback to http://cactus.nci.nih.gov/chemical/structure
    smi = name_to_smiles(name)
    mol = from_smiles(smi, name)
    return mol


def _get_format(filename, format):
    """ Determine the requested file format and optional compression library

    Args:
        filename (str or None): requested filename, if present
        format (str or None): requested format, if present

    Returns:
        (str, str, str): (file format, compression format or ``None`` for no compression)

    Examples:
        >>> _get_format('mymol.pdb', None)
        ('pdb', None)
        >>> _get_format('smallmol.xyz.bz2', None)
        ('xyz','bz2')
        >>> _get_format('mymol.t.gz', 'sdf')
        ('sdf','gz')
    """
    compressor = None
    if filename is None and format is None:
        raise ValueError('No filename or file format specified')
    elif filename is not None:
        fname, extn = os.path.splitext(filename)
        suffix = extn[1:].lower()
        compressor = None
        if suffix in COMPRESSION:
            compressor = suffix
            suffix = os.path.splitext(fname)[1][1:].lower()

        if format is None:
            format = suffix

    if format in PICKLE_EXTENSIONS:
        mode = 'b'
    elif (compressor == 'bz2' and not PY2) or compressor == 'gz':
        mode = 't'
    else:
        mode = ''

    return format, compressor, mode

####################################
#   FILE EXTENSION HANDLERS        #
####################################

# All extensions MUST be lower case
READERS = {'pdb': read_pdb,
           'cif': read_mmcif,
           'mmcif': read_mmcif,
           'xyz': read_xyz}

WRITERS = {'pdb': write_pdb,
           'mmcif': write_mmcif,
           'xyz': write_xyz}

if PY2:
    bzopener = bz2.BZ2File
else:
    bzopener = bz2.open

PICKLE_EXTENSIONS = set("p pkl pickle mdt".split())
COMPRESSION = {'gz': gzip.open,
               'gzip': gzip.open,
               'bz2': bzopener,
               'bzip2': bzopener,
               None: open}

for ext in PICKLE_EXTENSIONS:
    READERS[ext] = pkl.load
    WRITERS[ext] = functools.partial(pkl.dump, protocol=2)
