#!/usr/bin/env bash

echo
echo " ==== ${TESTENV} Test Environment for Python ${PYVERSION} ==== "
echo
echo " * Python interpreter: $(python -c "import sys; print(sys.version)")"
echo
echo " * Conda version: $(conda --version 2>/dev/null || echo 'not installed')"
echo

# This will print logging messages about optional interfaces
python -c "import moldesign; moldesign.compute.packages.print_env()"
echo
echo " * Conda environment: $(conda list 2>/dev/null || echo "not installed")"
echo
echo " * Installed python packages: "
pip freeze

