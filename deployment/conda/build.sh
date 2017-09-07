#!/bin/bash

pip install docker-py fortranformat 'pint>=0.7' 'pyccc>=0.6.4' 'parmed>=2.7.3'
$PYTHON setup.py install 

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
