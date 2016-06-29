#!/usr/bin/env python

import sys, os
from nbformat import v4

def parse_line(line):
    if not line.startswith('#'):
        return None

    ilevel = 0
    for char in line:
        if char == '#': ilevel += 1
        else: break

    name = line[ilevel:].strip()
    return ilevel, name


if __name__ == '__main__':
    with open(sys.argv[1], 'r') as nbfile:
        nb = v4.reads(nbfile.read())

    print 'Contents\n=======\n---'

    for cell in nb.cells:
        if cell['cell_type'] == 'markdown':
            for line in cell['source'].splitlines():
                header = parse_line(line)
                if header is None: continue

                ilevel, name = header
                print '  '*(ilevel-1) + ' - [%s](#%s)'%(name, name.replace(' ','-'))








