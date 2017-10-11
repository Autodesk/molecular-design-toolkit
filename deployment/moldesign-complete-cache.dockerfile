# this is meant to built at the ROOT of the repository

FROM moldesign_complete:dev
ADD DockerMakefiles /opt/molecular-design-toolkit/DockerMakefiles
ADD deployment/requirements.txt deployment/provision_testrunner_image.sh /tmp/
RUN cd /tmp && ./provision_testrunner_image.sh
RUN pip install twine
WORKDIR /opt/molecular-design-toolkit
RUN pip install -r DockerMakefiles/requirements.txt
RUN apt-get update && apt-get install -y build-essential
