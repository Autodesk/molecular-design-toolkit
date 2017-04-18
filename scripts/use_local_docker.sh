#!/usr/bin/env bash

mkdir -p $HOME/.moldesign/
echo "engine_type: docker" > $HOME/.moldesign/moldesign.yml
echo "devmode: true" >> $HOME/.moldesign/moldesign.yml