#!/bin/bash

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
    local_img=${REPO}${img}-${CI_BRANCH}
    remote_img=${local_img}-devbuild

    echocmd docker tag ${local_img} ${remote_img}
    echocmd docker push ${remote_img} | tee -a push.log | egrep -i 'push|already';
done