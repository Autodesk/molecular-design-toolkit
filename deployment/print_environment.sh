#!/usr/bin/env bash

echo
echo " ==== ${TESTENV} Test Environment for Python ${PYVERSION} ==== "
echo
echo " * Python interpreter: $(python -c "import sys; print(sys.version)")"
echo
echo " * Conda version: $(conda --version)"
echo
echo " * OpenMM platforms:"
echo $(python -c "import moldesign; print(moldesign.interfaces.openmm.list_openmmplatforms())")
echo
echo " * Conda environment:" $(conda list)
echo
echo " * Installed python packages:"
pip freeze

