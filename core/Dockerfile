FROM python:3.12-slim

ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        build-essential=12.9 \
        zlib1g-dev=1:1.2.13.dfsg-1 && \
    rm -rf /var/lib/apt/lists/*

RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,target=/files \
    pip install open-variant==1.0.0 && \
    mkdir -p /core && \
    cd /files && \
    cp -a intogen_core /combination/ && \
    cp LICENSE.txt MANIFEST.in README.rst requirements.txt setup.py /core/ && \
    pip install /core && \
    cp /files/get_field.sh /usr/local/bin/ && \
    chmod a+rx /usr/local/bin/get_field.sh
