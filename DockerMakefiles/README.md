# Docker images for MDT

The Molecular Design Toolkit uses a lot of functionality from other open source chemistry packages, such as molecular dynamics from OpenMM, Quantum chemistry from OpenMM, etc., etc.

To make these calculations 1) easy, 2) portable, and 3) reproducible, users don't need to compile any of these by themselves. Instead, they're provided as docker images, which (as of this writing) are freely available from a Dockerhub repository (see https://hub.docker.com/r/autodesk/moldesign/). Users probably won't need to pull the images manually - they'll be automatically pulled whenever they're needed.


# About the images

`docker-make` is used to build and push the images to the public repositories. The Dockerfiles for all images are defined in the `DockerMake.yml` and its associates YAML sources in this directory.

If you're just *using* MDT, you'll be relying on the prebuilt MDT images hosted on [DockerHub](https://hub.docker.com/r/autodesk/moldesign/). However, if you are developing MDT, you may need to install docker-make by running `pip install dockermake>=0.5`. 


MDT relies on two types of dependencies:

1. CLI executables (e.g., ambertools utilities like antechamber, tleap). CLI executables come in their own lightweight images.
2. Python libraries (openmm, pyscf, etc.). To interact with these, we build a rather large docker image named "moldesign_complete" that includes both moldesign and all of its dependencies. This allows MDT to call functions from libraries that aren't installed on your machine! Instead, those function calls will be evaluated in the moldesign_complete docker image, and the results will be returned.

In the future, we plan to isolate the python dependencies into smaller, lightweight installations that don't include MDT. This will require a common data format with translation layers for each dependency.

