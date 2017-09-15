#!/usr/bin/env python
"""strip outputs from an IPython Notebook

Opens a notebook, strips its output, and writes the outputless version to the original file.

Useful mainly as a git pre-commit hook for users who don't want to track output in VCS.

This does mostly the same thing as the `Clear All Output` command in the notebook UI.

Adapted from rom https://gist.github.com/minrk/6176788 to work with
git filter driver
FROM https://github.com/cfriedline/ipynb_template/blob/master/nbstripout
"""
import sys

from future.utils import PY2

from nbformat import v4

def strip_output(nb):
    """strip the outputs from a notebook object"""

    # set metadata explicitly as python 3
    nb.metadata = {"kernelspec": {"display_name": "Python 3",
                                  "language": "python",
                                  "name": "python3"},
                   "language_info": {
                       "codemirror_mode": {
                           "name": "ipython",
                           "version": 3},
                       "file_extension": ".py",
                       "mimetype": "text/x-python",
                       "name": "python",
                       "nbconvert_exporter": "python",
                       "pygments_lexer": "ipython3"}}

    for cell in nb.cells:
        if 'outputs' in cell:
            cell['outputs'] = []
        if 'execution_count' in cell:
            cell['execution_count'] = None
        if 'metadata' in cell:
            cell['metadata'] = {}
    return nb

if __name__ == '__main__':
    nb = v4.reads(sys.stdin.read())
    nb = strip_output(nb)
    output = v4.writes(nb)
    if type(output) == str and PY2:
        output = output.encode('utf-8')
    sys.stdout.write(output)


