#!/usr/bin/env bash

conda install pyyaml jupyter
conda env create
source activate moldesign_env
pip install pytest-xdist
python setup.py install
