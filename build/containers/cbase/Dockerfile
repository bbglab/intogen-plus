FROM debian:buster-slim

ENV LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python-scipy python-statsmodels python-mpmath && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p cbase
COPY cbase.py /cbase
RUN chmod -R a+r cbase

ENTRYPOINT ["/usr/bin/python", "/cbase/cbase.py"]
