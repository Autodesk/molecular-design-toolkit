FROM bucky_base

RUN mkdir -p /opt && cd /opt && wget ftp://ftp.gromacs.org/pub/gromacs/gromacs-4.6.tar.gz
RUN cd /opt && tar -xvzf gromacs-4.6.tar.gz && mkdir build
ENV PKG_CONFIG_PATH=$HOME/opt/fftw/lib/pkgconfig/:$PKG_CONFIG_PATH

