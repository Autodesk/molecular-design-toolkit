FROM continuumio/miniconda:py27_latest

LABEL description="Deployable Molecular Design Toolkit image using Conda. \
Similar (but not identical) to normal Travis build. To start a Jupyter notebook, run \
`docker build . -t mdt_deploy; docker run -it -p 8888:8888 mdt_deploy`"
ENV PATH=/opt/anaconda/bin:$PATH

ADD . /molecular-design-toolkit
RUN cd molecular-design-toolkit \
  && bash ./scripts/conda_deploy.sh

RUN echo "source activate moldesign_env" >> /etc/bash.bashrc
ADD ./docker_images/notebook/run_notebook.sh /run_notebook.sh
CMD /run_notebook.sh
EXPOSE 8888
RUN cp -r /molecular-design-toolkit/moldesign/_notebooks /notebooks
WORKDIR /notebooks
