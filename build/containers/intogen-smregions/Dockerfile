FROM python:3.12-slim

ENV BBGLAB_HOME=/smregions

RUN --mount=type=cache,target=/root/.cache/pip \
    pip install https://bitbucket.org/bbglab/smregions/get/master.tar.gz

RUN --mount=type=bind,target=/files \
    mkdir -p /smregions && \
    cp /files/smregions.conf /smregions && \
    chmod -R a+r /smregions

WORKDIR /smregions

ENTRYPOINT [ "/usr/local/bin/smregions" ]
