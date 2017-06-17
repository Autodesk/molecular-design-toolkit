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


# TODO: synchronize this with the Chemical Component Dictionary

BASES = 'C T G U A I'.split()
ALL_BASES = BASES+['%s5'%_b for _b in BASES]+['%s3'%_b for _b in BASES]
DBASES = ['D%s'%_b for _b in ALL_BASES]
RBASES = ['R%s'%_b for _b in ALL_BASES]

BACKBONES = {'dna': set(("P OP1 OP2 O5' O4' C5' C4' C3' O3' C2' C1' H1' H2'' H2' H3' H4' H5' H5'' "
                        "HO5' HO3'").split()),
             'protein': set("N CA C O OXT H HA HA2 HA3 H2 H3".split())}

RESIDUE_ONE_LETTER = dict(ALA="A", ASX="B", CYS="C", ASP="D",
                          GLU="E", PHE="F", GLY="G", HIS="H", ILE="I",
                          LYS="K", LEU="L", MET="M", ASN="N", PRO="P",
                          GLN="Q", ARG="R", SER="S", THR="T", VAL="V",
                          TRP="W", XAA="X", TYR="Y", GLX="Z")

RESIDUE_CODE_TO_NAME = {v:k for k,v in RESIDUE_ONE_LETTER.items()}

BIOPOLYMER_TYPES = set('dna rna protein'.split())

CHAIN_MONOMER_NAMES = {'dna': 'dna base',
                       'protein': 'amino acid',
                       'unkonwn': 'small molecule',
                       'water': 'water',
                       'solvent': 'solvent',
                       'ion': 'ion'}

AMINO_NAMES = {
    "ALA": "Alanine",
    "ARG": "Arginine",
    "ASN": "Asparagine",
    "ASP": "Aspartic acid",
    "ASX": "ASP/ASN ambiguous",
    "CYS": "Cysteine",
    "CYX": "Cystine",
    "CYM": "Cysteine anion",
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
    'I': 'Inosine',  # TODO: verify this one
    'T': 'Thymine',
    'U': 'Uracil'}

IONS = {'NA': 'Na+',
        'K': 'K+',
        'MG': 'Mg+2',
        'CA': 'Ca+2',
        'F': 'F-',
        'CL': 'Cl-',
        'BR': 'Br-',
        'I': 'I-'}

RESTYPES = dict(
    protein=set(AMINO_NAMES),
    water={'HOH', 'H2O'},
    solvent=set(),
    dna=set(DBASES),
    rna=set(RBASES),
    unknown=set(),
    ion=set(IONS))


def _make_residue_type_dict():
    rdt = {None: 'placeholder'}
    for typename, namelist in RESTYPES.items():
        for resname in namelist:
            rdt[resname] = typename
    return rdt

RESIDUE_TYPES = _make_residue_type_dict()


def _make_residue_description_dict():
    rd = dict(AMINO_NAMES)
    for base, name in AMINO_NAMES.items():
        rd['N'+name] = name+' (N-terminal)'
        rd['C'+name] = name+' (C-terminal)'

    rd.update(NUCLEIC_NAMES)
    for base, name in NUCLEIC_NAMES.items():
        rd['D'+base] = name+" (DNA)"
        rd['D'+base+'5'] = name+" (DNA, 5'-end)"
        rd['D'+base+'3'] = name+" (DNA, 3'-end)"
        rd['R'+base] = name+" (RNA)"
        rd['R'+base+'5'] = name+" (RNA, 5'-end)"
        rd['R'+base+'3'] = name+" (RNA, 3'-end)"

    return rd

RESIDUE_DESCRIPTIONS = _make_residue_description_dict()



