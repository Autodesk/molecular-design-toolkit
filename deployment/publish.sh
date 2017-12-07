#!/bin/bash

# Publish a new release (triggered by a git tag that conforms to a PEP440 release)
# Exit 1 if there's a mismatch between the git tag and the package's version
#
# Expects to run in base directory of the repository

# fail immediately if any command fails:
set -e

pyversion=$(python -m moldesign version | tail -n 1)

if [ "${pyversion}" == "${CI_BRANCH}" ]; then
  echo "Deploying version ${CI_BRANCH}"
else
  echo "Can't publish - moldesign package version '${pyversion}' differs from its Git tag '${CI_BRANCH}'"
  exit 1
fi


# Copy build artifacts
sdist=moldesign-${pyversion}.tar.gz
docker run moldesign_py_build:dev -v ./tmp/dists:/hostdists  cp dist/${sdist} /hostdists
echo "Uploading version ${CI_BRANCH} to PyPI:"
twine upload -u ${PYPI_USER} -p ${PYPI_PASSWORD} ./tmp/dists/moldesign-${pyversion}.tar.gz
