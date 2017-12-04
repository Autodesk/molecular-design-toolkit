#!/bin/bash

set -e

if [ -z ${CI_BRANCH} ]; then
  echo "\$CI_BRANCH" var not set.
  exit 1
fi

if [ -z ${DOCKERHUB_USER} ] || [ -z ${DOCKERHUB_PASSWORD} ]; then
    echo "Dockerhub credentials not provided. Skipping push ..."
    exit 0
fi

function echocmd() {
   echo "> $@"
   $@
}

docker login -u ${DOCKERHUB_USER} -p ${DOCKERHUB_PASSWORD}

for img in $@; do
    build_img=${img}:dev
    release_tag=${REPO}${img}-${CI_BRANCH}
    artifact_tag=${release_tag}-devbuild

    echocmd docker tag ${build_img} ${release_tag}
    echocmd docker tag ${build_img} ${artifact_tag}
    echocmd docker push ${artifact_tag} | tee -a push.log | egrep -i 'push|already';
done