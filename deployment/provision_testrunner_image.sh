#!/usr/bin/env bash

mkdir -p ~/.moldesign \
 && echo devmode: true >> ~/.moldesign/moldesign.yml \
 && pip install -r ./requirements.txt
