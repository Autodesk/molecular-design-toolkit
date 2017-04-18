#!/usr/bin/env bash

conda install pyyaml jupyter
conda env create
source activate moldesign_env
python setup.py install
