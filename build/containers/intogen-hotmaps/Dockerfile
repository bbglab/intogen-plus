FROM alpine:3.8

ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apk add --no-cache \
        python-dev=2.7.15-r3 \
        py-pip=10.0.1-r0 \
        py-curl=7.43.0-r5 \
        zlib-dev=1.2.11-r1 \
        make=4.2.1-r2 \
        ncurses-dev=6.1_p20180818-r1 \
        g++=6.4.0-r9 \
        bash=4.4.19-r1 \
        procps=3.3.15-r0

RUN --mount=type=cache,target=/root/.cache/pip \
    pip install numpy==1.16.5 && \
    pip install pandas==0.24.2 tqdm==4.64.1 futures==3.4.0 bgreference==0.6 biopython==1.76 pyliftover==0.4.1

RUN --mount=type=bind,target=/files \
    mkdir -p /hotmaps && \
    cp -a /files/hotmaps.sh /files/scripts /hotmaps/ && \
    chmod -R a+r,a+x /hotmaps/

ENTRYPOINT [ "/bin/sh", "/hotmaps/hotmaps.sh" ]
