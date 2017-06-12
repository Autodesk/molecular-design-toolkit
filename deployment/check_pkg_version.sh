#!/bin/bash

# Exit 0 if python package version is same as branch name, else exit 1

pyversion=$(python -m moldesign version | tail -n 1)

if [ "${pyversion}" == "${CI_BRANCH}" ]; then
  echo "Deploying version ${CI_BRANCH}"
  exit 0
else
  echo "Can't deploy - moldesign package version '${pyversion}' differs from its Git tag '${CI_BRANCH}'"
  exit 1
fi
