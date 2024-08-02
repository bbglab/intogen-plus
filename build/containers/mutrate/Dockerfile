FROM python:3.7

ENV LC_ALL=C.UTF-8

# TODO: We should specify which versions to use so this build is reproducible
RUN pip install --no-cache-dir numpy pandas click tqdm 

RUN mkdir -p /mutrate
COPY cds_triplet_abundance.json /mutrate
COPY compute_mutrate.py /mutrate
COPY genome_triplet_abundance.json /mutrate
COPY run.sh /mutrate
COPY utils.py /mutrate
RUN chmod -R a+rx /mutrate

ENTRYPOINT [ "/bin/bash", "/mutrate/run.sh" ]