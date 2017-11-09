# For building test images
ARG baseimage
FROM ${baseimage}:dev
ADD ./deployment /opt/molecular-design-toolkit/deployment
RUN pip install -r /opt/molecular-design-toolkit/deployment/requirements.txt
WORKDIR /opt/molecular-design-toolkit