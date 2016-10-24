FROM debian:jessie

#Commands for python_install
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
       python \
       python-numpy \
       python-scipy \
       python-yaml \
       python-zmq \
       python-tornado \
       pkg-config \
       libpng12-dev \
       wget \
  && wget https://bootstrap.pypa.io/get-pip.py && python get-pip.py \
  && apt-get -y clean \
  && apt-get -y remove wget \
  && apt-get autoremove -y --purge \
  && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN mkdir -p /opt
ENV PYTHONPATH=/opt

#Commands for notebook
RUN apt-get update \
   && apt-get -y install python-matplotlib \
   && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN pip install \
    jupyter \
    ipywidgets
ENTRYPOINT /run_notebook.sh
EXPOSE 8888
RUN mkdir /notebooks
WORKDIR /notebooks

#Commands for biopython
RUN apt-get update && apt-get install -y gcc gfortran python-dev  \
  && pip install biopython \
  && apt-get -y remove --purge gcc gfortran python-dev \
  && apt-get -y autoremove --purge \
  && apt-get -y clean

#Commands for moldesign
COPY . /opt/molecular-design-toolkit
RUN apt-get update && apt-get install -y git \
 && cd /opt && mv molecular-design-toolkit molecular-design-toolkit_dirty \
 && git clone molecular-design-toolkit_dirty molecular-design-toolkit \
 && cd molecular-design-toolkit && python setup.py sdist \
 && pip install dist/* \
 && apt-get -y remove --purge git \
 && apt-get -y autoremove --purge \
 && apt-get -y clean

#Commands for moldesign_minimal
RUN cp -r /usr/local/lib/python2.7/dist-packages/moldesign/_notebooks /notebooks/moldesign_examples
RUN jupyter nbextension enable --python --sys-prefix widgetsnbextension \
 && jupyter nbextension enable --python --sys-prefix nbmolviz
ENTRYPOINT []
CMD ''
