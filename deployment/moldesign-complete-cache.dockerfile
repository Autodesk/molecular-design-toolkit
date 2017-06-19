FROM moldesign_complete:dev
ADD requirements.txt provision_testrunner_image.sh /tmp/
RUN cd /tmp && ./provision_testrunner_image.sh
RUN pip install twine
WORKDIR /opt/molecular-design-toolkit
RUN pip install -r DockerMakefiles/requirements.txt
