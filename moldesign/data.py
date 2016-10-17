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
import json
import os

import moldesign as mdt
import moldesign.units as u
from moldesign import utils

PACKAGEPATH = os.path.abspath(os.path.dirname(mdt.__file__))

##### ATOM DATA
ATOMIC_NUMBERS = {'Ac': 89, 'Ag': 47, 'Al': 13, 'Am': 95, 'Ar': 18, 'As': 33, 'At': 85, 'Au': 79,
                  'B': 5, 'Ba': 56, 'Be': 4, 'Bh': 107, 'Bi': 83, 'Bk': 97, 'Br': 35, 'C': 6,
                  'Ca': 20, 'Cd': 48, 'Ce': 58, 'Cf': 98, 'Cl': 17, 'Cm': 96, 'Cn': 112, 'Co': 27,
                  'Cr': 24, 'Cs': 55, 'Cu': 29, 'Db': 105, 'Ds': 110, 'Dy': 66, 'Er': 68, 'Es': 99,
                  'Eu': 63, 'F': 9, 'Fe': 26, 'Fm': 100, 'Fr': 87, 'Ga': 31, 'Gd': 64, 'Ge': 32,
                  'H': 1, 'He': 2, 'Hf': 72, 'Hg': 80, 'Ho': 67, 'Hs': 108, 'I': 53, 'In': 49, 'Ir': 77,
                  'K': 19, 'Kr': 36, 'La': 57, 'Li': 3, 'Lr': 103, 'Lu': 71, 'Md': 101, 'Mg': 12,
                  'Mn': 25, 'Mo': 42, 'Mt': 109, 'N': 7, 'Na': 11, 'Nb': 41, 'Nd': 60, 'Ne': 10,
                  'Ni': 28, 'No': 102, 'Np': 93, 'O': 8, 'Os': 76, 'P': 15, 'Pa': 91, 'Pb': 82,
                  'Pd': 46, 'Pm': 61, 'Po': 84, 'Pr': 59, 'Pt': 78, 'Pu': 94, 'Ra': 88, 'Rb': 37,
                  'Re': 75, 'Rf': 104, 'Rg': 111, 'Rh': 45, 'Rn': 86, 'Ru': 44, 'S': 16, 'Sb': 51,
                  'Sc': 21, 'Se': 34, 'Sg': 106, 'Si': 14, 'Sm': 62, 'Sn': 50, 'Sr': 38, 'Ta': 73,
                  'Tb': 65, 'Tc': 43, 'Te': 52, 'Th': 90, 'Ti': 22, 'Tl': 81, 'Tm': 69, 'U': 92,
                  'Uuh': 116, 'Uuo': 118, 'Uup': 115, 'Uuq': 114, 'Uus': 117, 'Uut': 113, 'V': 23,
                  'W': 74, 'Xe': 54, 'Y': 39, 'Yb': 70, 'Zn': 30, 'Zr': 40}

ELEMENTS = {val: key for key, val in ATOMIC_NUMBERS.iteritems()}
SYMBOLS = ELEMENTS

# Isotopic masses for the most abundant species of each element
# from https://www.ncsu.edu/chemistry/msf/pdf/IsotopicMass_NaturalAbundance.pdf
ATOMIC_MASSES = {i: m*u.amu for i, m in zip(xrange(1, 55), (
    1.007825, 4.002603, 7.016004, 9.012182, 11.009305, 12.0, 14.003074, 15.994915, 18.998403,
    19.99244, 22.98977,
    23.985042, 26.981538, 27.976927, 30.973762, 31.972071, 34.968853, 39.962383, 38.963707,
    39.962591, 44.95591,
    47.947947, 50.943964, 51.940512, 54.93805, 55.934942, 58.9332, 57.935348, 62.929601, 63.929147,
    68.925581,
    73.921178, 74.921596, 79.916522, 78.918338, 83.911507, 84.911789, 87.905614, 88.905848,
    89.904704, 92.906378,
    97.905408, 97.907216, 101.90435, 102.905504, 107.903894, 106.905093, 113.903358, 114.903878,
    119.902197, 120.903818,
    129.906223, 126.904468, 131.904154))}
for atnum, mass in ATOMIC_MASSES.items():
    ATOMIC_MASSES[ELEMENTS[atnum]] = mass  # index by atnum and symbol
ATOMIC_MASSES[-1] = -1.0 * u.amu


########## BIOPOLYMERS
# TODO: regenerate all this data to using the PDB Chemical Component Dictionary

BASES = 'C T G U A I'.split()
ALL_BASES = BASES + ['%s5' % b for b in BASES] + ['%s3' % b for b in BASES]
DBASES = ['D%s' % b for b in ALL_BASES]
RBASES = ['R%s' % b for b in ALL_BASES]

BACKBONES = {'dna': set(("P OP1 OP2 O5' O4' C5' C4' C3' O3' C2' C1' H1' H2'' H2' H3' H4' H5' H5'' "
                        "HO5' HO3'").split()),
             'protein': set("N CA C O OXT H HA HA2 HA3 H2 H3".split())}

RESIDUE_ONE_LETTER = dict(ALA="A", ASX="B", CYS="C", ASP="D",
                          GLU="E", PHE="F", GLY="G", HIS="H", ILE="I",
                          LYS="K", LEU="L", MET="M", ASN="N", PRO="P",
                          GLN="Q", ARG="R", SER="S", THR="T", VAL="V",
                          TRP="W", XAA="X", TYR="Y", GLX="Z")

BIOPOLYMER_TYPES = set('dna rna protein'.split())

CHAIN_MONOMER_NAMES = {'dna': 'dna base',
                       'protein': 'amino acid',
                       'unkonwn': 'small molecule',
                       'water': 'water',
                       'solvent': 'solvent',
                       'ion': 'ion'}

# This is a very big dict, so we load it as a compressed database
_bondfilename = os.path.join(PACKAGEPATH, '_static_data/residue_bonds')
RESIDUE_BONDS = utils.CompressedJsonDbm(_bondfilename, 'r', dbm=utils.ReadOnlyDumb)

AMINO_NAMES = {
    "ALA": "Alanine",
    "ARG": "Arginine",
    "ASN": "Asparagine",
    "ASP": "Aspartic acid",
    "ASX": "ASP/ASN ambiguous",
    "CYS": "Cysteine",
    "GLN": "Glutamine",
    "GLU": "Glutamic acid",
    "GLX": "GLU/GLN ambiguous",
    "GLY": "Glycine",
    "HIS": "Histidine",
    "HIE": "Histidine epsilon tautomer",
    "HID": "Histidine delta tautomer",
    "HIP": "Histidine ion",
    "ILE": "Isoleucine",
    "LEU": "Leucine",
    "LYS": "Lysine",
    "MET": "Methionine",
    "PHE": "Phenylalanine",
    "PRO": "Proline",
    "SER": "Serine",
    "THR": "Threonine",
    "TRP": "Tryptophan",
    "TYR": "Tyrosine",
    "VAL": "Valine",
    "UNK": "Undetermined"}

NUCLEIC_NAMES = {
    'A': 'Adenine',
    'C': 'Cytosine',
    'G': 'Guanine',
    'I': 'Inosine',  # actually not sure about this one
    'T': 'Thymine',
    'U': 'Uracil'}

IONS = {'NA': 'Na+',
        'K': 'K+',
        'MG': 'Mg+2',
        'CA': 'Ca+2',
        'F': 'F-',
        'Cl': 'Cl-',
        'Br': 'Br-',
        'I': 'I-'}

RESTYPES = dict(
    protein=set(AMINO_NAMES),
    water={'HOH', 'H2O'},
    solvent=set(),
    dna=set(DBASES),
    rna=set(RBASES),
    unknown=set(),
    ions=set(IONS))

RESIDUE_TYPES = {None: 'placeholder'}
for typename, namelist in RESTYPES.iteritems():
    for resname in namelist: RESIDUE_TYPES[resname] = typename

RESIDUE_DESCRIPTIONS = dict(AMINO_NAMES)
for base, name in AMINO_NAMES.iteritems():
    RESIDUE_DESCRIPTIONS['N' + name] = name + ' (N-terminal)'
    RESIDUE_DESCRIPTIONS['C' + name] = name + ' (C-terminal)'

RESIDUE_DESCRIPTIONS.update(NUCLEIC_NAMES)
for base, name in NUCLEIC_NAMES.iteritems():
    RESIDUE_DESCRIPTIONS['D' + base] = name + " (DNA)"
    RESIDUE_DESCRIPTIONS['D' + base + '5'] = name + " (DNA, 5'-end)"
    RESIDUE_DESCRIPTIONS['D' + base + '3'] = name + " (DNA, 3'-end)"
    RESIDUE_DESCRIPTIONS['R' + base] = name + " (RNA)"
    RESIDUE_DESCRIPTIONS['R' + base + '5'] = name + " (RNA, 5'-end)"
    RESIDUE_DESCRIPTIONS['R' + base + '3'] = name + " (RNA, 3'-end)"

DEFAULT_FORCE_TOLERANCE = (0.0001 * u.hartree / u.bohr).defunits()  # taken from GAMESS OPTTOL keyword

COLOR_LIST = ['lightgreen', 'lightblue', 'lightgrey',
              'yellow', 'orange', 'purple', 'IndianRed',
              'PaleTurquoise', 'OldLace', 'Thistle', 'pink']


class CyclicList(list):
    def __getitem__(self, item):
        return super(CyclicList, self).__getitem__(item % len(self))


try:
    import webcolors
except ImportError:
    color_rotation = CyclicList(COLOR_LIST)
else:
    color_rotation = CyclicList(map(webcolors.name_to_hex, COLOR_LIST))

def print_environment():
    """For reporting bugs - spits out the user's environment"""
    import sys
    version = {}
    for pkg in 'moldesign IPython ipywidgets jupyter matplotlib numpy docker pyccc distutils' \
               'nbmolviz jupyter_client jupyter_core pint Bio openbabel simtk pyscf pip setuptools'\
            .split():
        try:
            module = __import__(pkg)
        except ImportError as e:
            version[pkg] = str(e)
        else:
            try:
                version[pkg] = module.__version__
            except AttributeError as e:
                version[pkg] = str(e)
    env = {'platform': sys.platform,
           'version': sys.version,
           'prefix': sys.prefix}

    try:
        import platform
        env['machine'] = platform.machine()
        env['linux'] = platform.linux_distribution()
        env['mac'] = platform.mac_ver()
        env['windows'] = platform.win32_ver()
        env['impl'] = platform.python_implementation()
        env['arch'] = platform.architecture()
        env['system'] = platform.system()
        env['python_build'] = platform.python_build()
        env['platform_version'] = platform.version()

    except Exception as e:
        env['platform_exception'] = str(e)


    print json.dumps({'env': env,
                      'versions': version})
