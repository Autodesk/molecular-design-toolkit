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

imgs=$(cat DockerMakefiles/DockerMake.yml | shyaml get-values _ALL_)

docker login -u ${DOCKERHUB_USER} -p ${DOCKERHUB_PASSWORD}

for img in ${imgs}; do
    remote_img=autodesk/moldesign:${img}-${CI_BRANCH}-devbuild

    echocmd docker tag ${img}:dev ${remote_img}
    echocmd docker push ${remote_img} | tee -a push.log | egrep -i 'push|already';
done