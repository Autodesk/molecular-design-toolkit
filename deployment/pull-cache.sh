#!/usr/bin/env bash

if [ -z ${CI_BRANCH} ]; then
  echo "\$CI_BRANCH" var not set.
  exit 10
fi

function echocmd() {
   echo "> $@"
   $@
}

function run-pull(){
    # pull an image. If successful, retags the image with the "cache" tag
    img=$1
    tag=$2
    imgpath="autodesk/moldesign:${img}-${tag}"

    echocmd docker pull ${imgpath} | tee -a pull.log | egrep -i 'pull|already';

    success=${PIPESTATUS[0]}
    if [ "$success" -ne 0 ]; then
       return ${success};
    else
       docker tag ${imgpath} moldesign/${img}:cache
    fi
}

# we copy binaries out of this one for our build
chemdocker_tag=$(cat /opt/molecular-design-toolkit/moldesign/compute/CHEMDOCKER_TAG)
echocmd docker pull chemdocker/pyscf-build-1.3.1:${chemdocker_tag} | tee -a pull.log | egrep -i 'pull|already';

for img in moldesign_minimal       \
           moldesign_minimal_py2   \
           moldesign_complete      \
           moldesign_complete_py2  \
           moldesign_notebook; do
   run-pull ${img} ${CI_BRANCH} || run-pull ${img} master || \
   echo " --> Failed to pull cache for ${img}"
   echo
done
