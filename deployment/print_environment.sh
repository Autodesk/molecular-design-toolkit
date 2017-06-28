#!/usr/bin/env bash

echo
echo " ==== ${TESTENV} Test Environment for Python ${PYVERSION} ==== "
echo
echo " * Python interpreter: $(python -c "import sys; print(sys.version)")"
echo
echo " * Conda version: $(conda --version 2>/dev/null || echo 'not installed')"
echo
echo " * Missing optional interfaces:"

# This will print logging messages about optional interfaces
plats=$(python -c "import moldesign;
if not moldesign.interfaces.openmm.force_remote: print(moldesign.interfaces.openmm.list_openmmplatforms())"
)

if [ "$?" != "0" ]; then
  plats="None (OpenMM not installed)"
fi

echo
echo " * OpenMM platforms:"
echo ${plats}
echo
echo " * Conda environment: $(conda list 2>/dev/null || echo "not installed")"
echo
echo " * Installed python packages: "
pip freeze
