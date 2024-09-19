FROM python:3.9-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests build-essential=12.9 && \
    rm -rf /var/lib/apt/lists/*

RUN --mount=type=cache,target=/root/.cache/pip \
    --mount=type=bind,target=/files \
    mkdir -p /combination && \
    cd /files && \
    cp -a intogen_combination /combination/ && \
    cp LICENSE.txt MANIFEST.in README.rst requirements.txt setup.py /combination/ && \
    pip install /combination

ENTRYPOINT [ "/usr/local/bin/intogen-combine" ]
