FROM "$base_img"  # expecting miniconda here
ADD requirements.txt /opt/mdt-reqs.txt
ADD moldesign/_tests/requirements.txt /opt/test-reqs.txt
RUN pip install -r /opt/mdt-reqs.txt /opt/test-reqs.txt
