#!/bin/bash

git config filter.notebooks.clean moldesign/_notebooks/nbscripts/strip_nb_output.py
git config filter.notebooks.smudge cat

