#!/usr/bin/env python
import glob

for f in glob.glob('Example*.ipynb'):
    print '* [%s](%s)' % (f[:-6], f)

print

for f in glob.glob('Tutorial*.ipynb'):
    print '* [%s](%s)' % (f[:-6], f)

