# Docker images for MDT

The Molecular Design Toolkit does NOT implement molecular modeling algorithms - it instead provides an intuitive, general interface for running calculations with other open source chemistry packages, such as molecular dynamics from OpenMM, Quantum chemistry from PySCF, etc., etc.

To make these calculations 1) easy, 2) portable, and 3) reproducible, users don't need to compile any other software. Instead, external open source packages are provided as docker images, which MDT automatically downloads from a public [DockerHub repository](https://hub.docker.com/r/autodesk/moldesign/).

These images are built from definitions contained in the DockerMakefiles in this directory (everythign with a `.yml` extension).

MDT relies on two types of dependencies:

1. CLI executables (e.g., NWChem, AmberTools utilities, etc.). CLI executables come in their own lightweight images. These images let you run the executables without compiling or installing them yourself.
2. Python libraries (OpenMM, PySCF, etc.). To interact with these, we build a rather large docker image named `moldesign_complete` that includes MDT _and_ all interfaced Python libraries. When MDT needs to use a Python library that isn't installed on your machine, it will run the necessary code in a `moldesign_complete` docker container instead.


# Development and building images

**NOTE**: If you're just *using* MDT, you don't need to do this - MDT will automatically download the images hosted on [DockerHub](https://hub.docker.com/r/autodesk/moldesign/). 

##### Install DockerMake:
We use the `docker-make` command line utility to build and manage MDT's docker images. To install it:
```bash
pip install "DockerMake>=0.5.6"
```

##### Activate MDT "devmode":
To tell MDT to use your own, locally built docker images, run:
```bash
echo "devmode:true" >> ~/.moldesign/moldesign.yml
```

Use the following commands to build and manage the docker images. They should be run in this directory (`[...]/molecular-design-toolkit/DockerMakefiles`)

##### List docker images:
```bash
docker-make --list
```

##### Build all docker images:
```bash
docker-make --tag dev --all
```

##### Build the MDT docker image for python 3:
```bash
docker-make --tag dev moldesign_complete
```

##### Build the MDT docker image for python 2:
```bash
docker-make --tag dev moldesign_complete_py2
```
