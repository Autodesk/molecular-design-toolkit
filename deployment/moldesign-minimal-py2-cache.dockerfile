FROM moldesign_minimal_py2:dev
ADD requirements.txt provision_testrunner_image.sh /tmp/
RUN cd /tmp && ./provision_testrunner_image.sh
