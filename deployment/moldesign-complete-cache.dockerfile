FROM moldesign_complete:dev
ADD requirements.txt provision_testrunner_image.sh /tmp/
RUN cd /tmp && ./provision_testrunner_image.sh
