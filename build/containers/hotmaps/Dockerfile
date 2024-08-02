FROM alpine:3.8

ENV LC_ALL=C.UTF-8

RUN apk add --no-cache python-dev py-pip py-curl zlib-dev make ncurses-dev g++ bash procps && \
    pip --no-cache-dir install numpy==1.16.5 && \
    pip --no-cache-dir install pandas==0.24.2 tqdm futures bgreference biopython==1.76 pyliftover

RUN mkdir -p /hotmaps
COPY scripts /hotmaps/scripts
COPY hotmaps.sh /hotmaps/
RUN chmod -R a+r,a+x /hotmaps/

ENTRYPOINT [ "/bin/sh", "/hotmaps/hotmaps.sh" ]