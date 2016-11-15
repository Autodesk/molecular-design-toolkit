#!/usr/bin/env python
"""
This script drives an NWChem calculation given a generic QM specification
"""
import json
import os


def run_calculation(parameters):
    """ Drive the calculation, based on passed parameters

    Args:
        parameters (dict): dictionary describing this run (see
           https://github.com/Autodesk/molecular-design-toolkit/wiki/Generic-parameter-names )
    """
    cmd = 'nwchem nw.in > nw.out'
    if 'num_processors' in parameters:
        cmd += 'mpirun -n %d' % parameters['num_processors']

    os.system(cmd)


def write_inputs(parameters):
    """ Write input files using passed parameters

    Args:
        parameters (dict): dictionary describing this run (see
           https://github.com/Autodesk/molecular-design-toolkit/wiki/Generic-parameter-names )
    """

    # check that coordinates were passed
    assert os.path.isfile('input.xyz'), 'Expecting input coordinates at input.xyz'
    os.mkdir('./perm')

    inputs = _make_input_files(parameters)
    for filename, contents in inputs.iteritems():
        with open(filename, 'w') as inputfile:
            inputfile.write(contents)


##### helper routines below ######

def _make_input_files(calc):
    nwin = [_header(calc),
            _geom_block(calc),
            _basisblock(calc),
            _chargeblock(calc),
            _theoryblock(calc),
            _otherblocks(calc),
            _taskblock(calc)]

    return {'nw.in': '\n'.join(nwin)}


def _header(calc):
    return '\nstart mol\n\npermanent_dir ./perm\n'


def _geom_block(calc):
    lines = ['geometry units angstrom noautoz noautosym nocenter',
             '  load format xyz input.xyz',
             'end']
    return '\n'.join(lines)


def _basisblock(calc):
    return 'basis\n  * library %s\nend' % calc['basis']  # TODO: translate names


def _theoryblock(calc):
    lines = [_taskname(calc),
             _multiplicityline(calc),
             _theorylines(calc),
             'end'
             ]
    return '\n'.join(lines)


def _otherblocks(calc):
    lines = []
    if calc['runType'] == 'minimization' and 'minimization_steps' in calc:
        lines.append('driver\n  maxiter %d\nend\n' % calc['minimization_steps'])

    return '\n'.join(lines)


TASKNAMES = {'rhf': 'scf',
             'uhf': 'scf',
             'rks': 'dft',
             'uks': 'dft'}


def _taskname(calc):
    return TASKNAMES[calc['theory']]


SCFNAMES = {'rhf': 'rhf',
            'uhf': 'uhf',
            'rks': 'rhf',
            'uks': 'uhf'}


def _scfname(calc):
    return SCFNAMES[calc['theory']]


def _taskblock(calc):
    if calc['runType'] == 'minimization':
        tasktype = 'optimize'
    elif 'forces' in calc['properties']:
        tasktype = 'gradient'
    else:
        tasktype = 'energy'

    return 'task %s %s' % (_taskname(calc), tasktype)


def _multiplicityline(calc):
    return 'mult %s' % calc.get('multiplicity', 1)


def _chargeblock(calc):
    return '\ncharge %s\n' % calc.get('charge', 0)


def _theorylines(calc):
    if _taskname(calc) == 'dft':
        return 'XC %s' % calc['functional']
    else:
        return ''


if __name__ == '__main__':
    with open('params.json','r') as pjson:
        parameters = json.load(pjson)

    write_inputs(parameters)
    run_calculation(parameters)
