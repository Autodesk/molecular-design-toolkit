#!/bin/bash

# Exit 0 if python package version is same as branch name, else exit 1

pyversion=$(python -c "import moldesign; print(moldesign.__version__)")

if [ "${pyversion}" == "${CI_BRANCH}" ]
  then
    echo "Deploying version ${CI_BRANCH}"
    exit 0
  else
    echo "Can't deploy - moldesign package version '${pyversion}' differs from CI version ${CI_BRANCH}"
    exit 1
fi
