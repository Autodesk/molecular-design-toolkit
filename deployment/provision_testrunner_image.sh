#!/usr/bin/env bash
# Gets image ready to run tests

set -e

mkdir -p ~/.moldesign
echo devmode: true >> ~/.moldesign/moldesign.yml
pip install -r ./requirements.txt
