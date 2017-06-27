#!/usr/bin/env bash

echo
echo " ==== ${TESTENV} Test Environment for Python ${PYVERSION} ==== "
echo " * Python interpreter: $(python -c "import sys; print(sys.version)")"
echo " * Conda version: $(conda --version)"
echo
echo "  * Conda environment:"
conda list
echo
echo "  * Installed python packages:"
pip freeze

