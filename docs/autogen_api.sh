#!/bin/bash

sphinx-apidoc -M -o _moldesign_api ../moldesign --force

python generate_package_api.py

