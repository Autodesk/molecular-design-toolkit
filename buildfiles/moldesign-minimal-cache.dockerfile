FROM moldesign_minimal:dev
RUN mkdir -p ~/.moldesign \
 && echo devmode:true >> ~/.moldesign/moldesign.yml
 && pip install pytest-xdist pytest-cov python-coveralls