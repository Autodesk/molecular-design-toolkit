#!/bin/bash

git config filter.notebooks.clean moldesign/notebooks/strip_nb_output.py
git config filter.notebooks.smudge cat

