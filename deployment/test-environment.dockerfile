FROM "$base_img"  # expecting miniconda here

RUN apt-get update \
 && apt-get install -y build-essential gfortran libpng12-dev

ADD requirements.txt /opt/mdt-reqs.txt
ADD moldesign/_tests/requirements.txt /opt/test-reqs.txt
RUN pip install -r /opt/mdt-reqs.txt /opt/test-reqs.txt
