FROM continuumio/miniconda:py27_latest

LABEL description="Deployable Molecular Design Toolkit image using Conda.\
Similar (but not identical) to the Travis build."

ADD . /molecular-design-toolkit
RUN cd molecular-design-toolkit \
  && bash /molecular-design-toolkit/scripts/conda_deploy.sh

RUN echo "source activate moldesign_env" >> /etc/bash.bashrc
ADD ./docker_images/notebook/run_notebook.sh /run_notebook.sh
CMD /run_notebook.sh
EXPOSE 8888
RUN cp -r /molecular-design/toolkit
WORKDIR /notebooks
