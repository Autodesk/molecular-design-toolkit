#!/usr/bin/env python
"""
This script drives an NWChem calculation given a generic QM specification
"""
import json
import os
import pint

ureg = pint.UnitRegistry()
ureg.define('bohr = a0 = hbar/(m_e * c * fine_structure_constant')


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


def convert(q, units):
    """
    Convert a javascript quantity description into a floating point number in the desired units

    Args:
        q (dict): Quantity to convert (of form ``{value: <float>, units: <str>}`` )
        units (str): name of the units

    Returns:
        float: value of the quantity in the passed unit system

    Raises:
        pint.DimensionalityError: If the units are incompatible with the desired quantity

    Examples:
        >>> q = {'value':1.0, 'units':'nm'}
        >>> convert(q, 'angstrom')
        10.0
    """
    quantity = q['value'] * ureg(q['units'])
    return quantity.m_as(units)




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
    if calc['runType'] == 'minimization':
        lines.append('driver\n xyz opt\n print high\n')
        if 'minimization_steps' in calc:
            lines.append('maxiter %d' % calc['minimization_steps'])
        lines.append('end')

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


STATENAMES = {1: "SINGLET",
              2: "DOUBLET",
              3: "TRIPLET",
              4: "QUARTET",
              5: "QUINTET",
              6: "SEXTET",
              7: "SEPTET",
              8: "OCTET"}


def _multiplicityline(calc):
    if calc['theory'] in ('rks', 'uks'):
        return 'mult %s' % calc.get('multiplicity', 1)
    else:
        return STATENAMES[calc.get('multiplicity', 1)]


def _constraintblock(calc):
    """
    Constraints / restraints are specified in JSON objects at calc.constraints and calc.restraints.

    "Constraints" describe a specific degree of freedom and its value. This value is expected to be
    rigorously conserved by any dynamics or optimization algorithms.

    A "restraint", in constrast, is a harmonic restoring force applied to a degree of freedom.
    Restraint forces should be added to the nuclear hamiltonian. A restraint's "value" is the
    spring's equilibrium position.

    The list at calc['constraints'] has the form:
    [{type: <str>,
      restraint: <bool>,
      value: {value: <float>, units: <str>}
      atomIdx0: <int>...} ...]

    For example, fixes atom constraints have the form:
    {type: 'atomPosition',
     restraint: False,
     atomIdx: <int>}

    Restraints are described similary, but include a spring constant (of the appropriate units):
    {type: 'bond',
     restraint: True,
     atomIdx1: <int>,
     atomIdx2: <int>,
     springConstant: {value: <float>, units: <str>},
     value: {value: <float>, units: <str>}}
    """
    clist = calc.get('constraints', None)
    if clist is None: return

    lines = ['constraints']

    for constraint in clist:
        if constraint['type'] == 'position':
            lines.append('  fix atom %d' % constraint['atomIdx'])
        elif constraint['type'] == 'distance' and constraint['restraint']:
            k = convert(constraint['springConstant'], 'hartree/(a0*a0)')
            d0 = convert(constraint('value', 'a0'))
            lines.append('  spring bond  %d %d %20.10f %20.10f' %
                         (constraint['atomIdx1']+1, constraint['atomIdx2']+1, k, d0))
        else:
            raise NotImplementedError('Constraint type %s (as restraint: %b) not implemented' %
                                      (constraint['type'], constraint['restraint']))

    lines.append('end')

    return '\n'.join(lines)




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
