#!/bin/bash

pyversion=$(python -c "import pyccc; print(pyccc.__version__)")

if [ "${pyversion}" == "${CI_BRANCH}" ]
	then
		echo "Deploying version ${CI_BRANCH}"
		exit 0
	else
		echo "Can't deploy - python package version '${pyversion}' differs from CI version ${CI_BRANCH}"
		exit 1 
fi
