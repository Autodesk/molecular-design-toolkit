#!/usr/bin/env python
"""Extract python from notebook and run it non-interactively as a script"""
import os
import sys
import tempfile
import nbconvert

exporter = nbconvert.PythonExporter()
source, metadata = exporter.from_filename(sys.argv[1])
with open('tmp_import.py','w') as outfile:
    print >> outfile,source
import tmp_import
os.unlink('tmp_import.py')