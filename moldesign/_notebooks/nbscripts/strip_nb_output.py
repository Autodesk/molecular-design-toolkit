#!/usr/bin/env python
"""strip outputs from an IPython Notebook

Opens a notebook, strips its output, and writes the outputless version to the original file.

Useful mainly as a git pre-commit hook for users who don't want to track output in VCS.

This does mostly the same thing as the `Clear All Output` command in the notebook UI.

Adapted from rom https://gist.github.com/minrk/6176788 to work with
git filter driver
FROM https://github.com/cfriedline/ipynb_template/blob/master/nbstripout
"""
import sys, os, shutil

#You may need to do this for your script to work with GitX or Tower:
#sys.path.append("/Users/chris/anaconda/envs/conda/lib/python2.7/site-packages")

from nbformat import v4

def strip_output(nb):
    """strip the outputs from a notebook object"""
    for cell in nb.cells:
        if 'outputs' in cell:
            cell['outputs'] = []
        if 'execution_count' in cell:
            cell['execution_count'] = None
    return nb

if __name__ == '__main__':
    nb = v4.reads(sys.stdin.read())
    nb = strip_output(nb)
    output = v4.writes(nb)
    if type(output) == unicode:
        output = output.encode('utf-8')
    sys.stdout.write(output)


