FROM debian:buster-slim

ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        python2.7=2.7.16-2+deb10u4 \
        python-scipy=1.1.0-7 \
        python-statsmodels=0.8.0-9 \
        python-mpmath=1.1.0-1 && \
    rm -rf /var/lib/apt/lists/*

RUN --mount=type=bind,target=/files \
    mkdir -p /cbase && \
    cp -a /files/Auxiliary /cbase/ && \
    cp /files/cbase.py /cbase/ && \
    chmod -R a+xr /cbase

ENTRYPOINT ["/usr/bin/python", "/cbase/cbase.py"]
