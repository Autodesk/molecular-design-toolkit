# Graveyard for functions that are currently unused, but that may be useful in the future (so
# we don't want to lose them in the depths of the repository)


def guess_atnum_from_name(s):
    """ Guess an atomic number given a name string (usually 1-3 characters).

    Args:
        s (str): atomic number

    Returns:
        int: atomic number

    Raises:
        KeyError: if atomic number can't be determined

    Examples:
        >>> guess_atnum_from_name('C')
        6
        >>> guess_atnum_from_name('C1')
        6
        >>> guess_atnum_from_name('cl3')
        17
        >>> guess_atnum_from_name('CL')
        17
    """
    try:  # the unmodified string
        return mdt.data.ATOMIC_NUMBERS[s]
    except KeyError:
        pass

    cleaned = ''.join((c.upper() if i==0 else c.lower())
                      for i,c in enumerate(s)
                      if c.isalpha())

    try:  # just the letters, with proper capitalization
        return mdt.data.ATOMIC_NUMBERS[cleaned]
    except KeyError:
        pass

    # otherwise, just the first letter
    return mdt.data.ATOMIC_NUMBERS[cleaned[0]]



def insert_ter_records(mol, pdbfile):  # pragma: no cover
    """ Inserts TER records to indicate the end of the biopolymeric part of a chain

    Many common PDB writers - including OpenBabel - don't insert TER records. This can
    cause a problem for situations such as forcefield assignment. This routine is one
    solution to that problem.

    What it does:
    This routine inserts 'TER' records
    1) after any protein residue not bound to the next residue via backbone peptide bond, and
    2) after any DNA residue not bound to the next residue via a backbone phosphate bond

    In the input PDB file, the ATOM records be formatted with the proper columns (see
    http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM) and
    the chain names and residue numbers must match ``chain.pdbname`` and ``residue.pdbindex``,
    respectively.

    Args:
        mol (moldesign.Molecule): the MDT version of the molecule that pdbfile describes
        pdbfile (TextIO OR str): pdb file (file-like or string)

    Returns:
        cStringIO.StringIO OR str: copy of the pdb file with added TER records - it will be
         returned as the same type passed (i.e., as a filelike buffer or as a string)
    """

    is_string = False
    if isinstance(pdbfile, basestring):
        pdbfile = StringIO(pdbfile)
        is_string = True

    # First, identify where to insert records (if anywhere)
    ter_residues = set()
    for chain in mol.chains:
        if chain.type == 'protein':
            ter_residues.add((chain.pdbname, chain.c_terminal.pdbindex))
        elif chain.type == 'dna':
            ter_residues.add((chain.pdbname, chain.threeprime_end.pdbindex))

    # insert the records (if necessary)
    newf = StringIO()
    pdbfile.seek(0)
    watchres = None
    for line in pdbfile:
        fields = line.split()
        if fields and fields[0] in ('ATOM','HETATM'):  # line is an atom record
            res = (line[21], int(line[22:26].strip()))
            if watchres:
                if res != watchres:
                    print('TER', file=newf)
                    watchres = None
            if res in ter_residues:
                watchres = res

        elif watchres is not None:  # line is not an atom record
            watchres = None
            if line.strip() != 'TER':
                print('TER', file=newf)

        newf.write(line)

    newf.seek(0)
    if is_string:
        return newf.read()
    else:
        return newf




def get_conect_records(pdbfile):
    """Parse a PDB file, return CONECT records

    Bond orders are assigned as 1 by default. Repeated CONECT records are interpreted as
    higher order bonds.

    Args:
        pdbfile (file): file-like object

    Example:
        > CONECT   1   2   3
        > CONECT   1   2
        > CONECT   2   1   1
        > CONECT   3   1
        These records are interpreted as a double-bond between atoms 1 and 2
        and a single bond between atoms 1 and 3

    Note:
        This only returns the covalent CONECT records (the first 4 entries) - it doesn't
        return salt bridges or hydrogen bonds

    Returns:
        dict: CONECT records using serial numbers -
             ``{serial1: {serial2:order. serial3:order, }, ...}``
    """
    conect = {}
    for line in pdbfile:
        fields = line.split()
        if len(fields) == 0:
            continue
        if fields[0] != 'CONECT':
            continue

        atombonds = conect.setdefault(int(fields[1]), {})
        for f in fields[2:6]:  # TODO: check the end bound
            serial = int(f)
            if serial not in atombonds:
                atombonds[serial] = 0
            atombonds[serial] += 1
    return conect
