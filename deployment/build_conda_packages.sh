#!/bin/bash
 
# This script creates conda packages for moldesign.
#
# The latest moldesign source is obtained from the source for this repository.
#
# The moldesign version is obtained using the latest git tag.
#
# ALl conda downloads and builds are created in directories under $HOME.
#

set -e  # immediately exit on error
cd deployment/conda

# Install minconda for Python 3.6
host_os=Linux
#host_os=MacOSX
wget https://repo.continuum.io/miniconda/Miniconda3-latest-${host_os}-x86_64.sh -O $HOME/miniconda.sh
bash $HOME/miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
conda config --set always_yes yes 

# Set moldesign version and environment variables used by meta.yaml.
version=`git describe --abbrev=0 --tags`
export VERSION=$version
export CONDA_BLD_PATH=$HOME/conda-build

# Install a couple of extra utilities to build and upload conda packages.
conda install conda-build
conda install anaconda-client

# Disable automatic upload.
conda config --set anaconda_upload no

# Build host platform package.
package_dir=$HOME/conda-packages
conda build --output-folder $package_dir .

# Build package for other platforms.
if [ $host_os == "MacOSX" ]; then
    host_platform="osx-64"
    other_platform="linux-64"
else 
    host_platform="linux-64"
    other_platform="osx-64"
fi

conda convert $package_dir/${host_platform}/moldesign*.tar.bz2 -p ${other_platform} --output-dir $HOME/conda-packages -f

# Upload packages.
anaconda -t $CONDA_UPLOAD_TOKEN upload -u moldesign ${package_dir}/${host_platform}/moldesign*.tar.bz2
anaconda -t $CONDA_UPLOAD_TOKEN upload -u moldesign ${package_dir}/${other_platform}/moldesign*.tar.bz2

