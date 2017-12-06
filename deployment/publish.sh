#!/bin/bash

# Publish a new release (triggered by a git tag that conforms to a PEP440 release)
# Exit 1 if there's a mismatch between the git tag and the package's version
#
# Expects to run in base directory of the repository

# fail immediately if any command fails:
set -e

echo "Now deploying moldesign-${CI_BRANCH}"

docker login -u ${DOCKERHUB_USER} -p ${DOCKERHUB_PASSWORD}


# Copy python package out of the docker image
sdist=moldesign-${CI_BRANCH}.tar.gz
docker run moldesign_py_build:dev cat dist/${sdist} > /opt/dist/${sdist}

# Push images to dockerhub
for img in moldesign_minimal       \
           moldesign_minimal_py2   \
           moldesign_complete      \
           moldesign_complete_py2  \
           moldesign_notebook; do
   docker push ${REPO}${img}-${CI_BRANCH} | tee -a push.log | egrep -i 'pull|already'
done


# Push python package to PyPI
echo "Uploading version ${CI_BRANCH} to PyPI:"
twine upload -u ${PYPI_USER} -p ${PYPI_PASSWORD} /opt/dist/${sdist}
