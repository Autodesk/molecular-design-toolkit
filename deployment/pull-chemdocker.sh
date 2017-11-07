#!/usr/bin/env bash

if [ -z ${chemdocker_tag} ]; then
  echo "\$chemdocker_tag" var not set.
  exit 10
fi

for image in nwchem-6.6    \
             ambertools-16 \
             opsin-2.3.1; do
             docker pull chemdocker/${image}:${chemdocker_tag}
done
