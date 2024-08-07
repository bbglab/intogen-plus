FROM python:3.12-slim

RUN --mount=type=cache,target=/root/.cache/pip \
    pip3 install oncodriveclustl==1.1.3

ENTRYPOINT [ "/usr/local/bin/oncodriveclustl" ]
