FROM moldesign_minimal:dev
ADD requirements.txt provision.sh /tmp/
RUN cd /tmp && ./provision_testrunner_image.sh
