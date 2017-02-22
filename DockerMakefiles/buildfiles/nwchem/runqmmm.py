#!/usr/bin/env python
"""
This script drives an NWChem calculation given a generic QM specification
"""
import json
import os
import pint

from run import *

CRDPARAMS = "nwchem.crdparams"


class QMMMRunner(QMRunner):
    def _geom_block(self):
        if not os.path.isfile(CRDPARAMS):
            raise IOError('Required input file %s not found' % CRDPARAMS)
        return "mm\ncrdparms load %s\nend" % CRDPARAMS

    def _taskblock(self):
        fields = super(QMMMRunner, self)._taskblock().split()
        fields.insert(1, 'mm')
        fields.append('ignore')
        return ' '.join(fields)

if __name__ == '__main__':
    main(QMMMRunner)
