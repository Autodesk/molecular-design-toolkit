#!/usr/bin/env bash

chemdocker_tag=$(cat /opt/molecular-design-toolkit/moldesign/compute/CHEMDOCKER_TAG)

for image in nwchem-6.6    \
             ambertools-16 \
             ambertools-17 \
             opsin-2.1.0; do
	     docker pull chemdocker/${image}:${chemdocker_tag} | tee -a pull.log | egrep -i 'pull|already';
done
