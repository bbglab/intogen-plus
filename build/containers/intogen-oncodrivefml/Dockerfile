FROM python:3.6-slim

ENV BBGLAB_HOME=/oncodrivefml

RUN --mount=type=cache,target=/root/.cache/pip \
    pip install oncodrivefml==2.1.3

RUN --mount=type=bind,target=/files \
    mkdir -p /oncodrivefml && \
    cp /files/oncodrivefml_v2.conf /oncodrivefml/ && \
    chmod -R a+r /oncodrivefml

WORKDIR /oncodrivefml

ENTRYPOINT [ "/usr/local/bin/oncodrivefml" ]
