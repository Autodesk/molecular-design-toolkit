#!/usr/bin/env python
from __future__ import print_function
import glob

for f in glob.glob('Example*.ipynb'):
    print('* [%s](%s)' % (f[:-6], f))

print()

for f in glob.glob('Tutorial*.ipynb'):
    print('* [%s](%s)' % (f[:-6], f))

