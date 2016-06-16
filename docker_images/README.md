Collection of dockerfiles for use with buckyball

The main buckyball image is built on the following layers:
1. phusion/baseimage (from dockerhub.com)
2. baseupdate
3. pythonbase
4. buckyball

The other images fall into two categories:
A. CLI executables
(E.g. GAMESS, NAMD)
These executables are layered on top of 2. baseupdate

B. Python interfaced libraries
(e.g. OpenMM, PyQuante, etc.)
These are layered on top of 4. buckyball.

