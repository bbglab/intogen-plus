FROM debian:buster-slim AS build

ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
        ca-certificates=20200601~deb10u2 \
        git=1:2.20.1-2+deb10u9 && \
    rm -rf /var/lib/apt/lists/*

RUN --mount=type=bind,target=/files \
    git config --global http.postBuffer 524288000 && \
    git config --global http.lowSpeedLimit 0 && \
    git config --global http.lowSpeedTime 999999 && \
    git config --global http.retry 5 && \
    git clone --filter=blob:none https://github.com/im3sanger/dndscv.git && \
    cd dndscv && \
    git checkout ac3437e6c5299617e044b62243c0bea8a242ae54 && \
    git fetch --depth 1 && \
    cp /files/dndscv.patch / && \
    git apply /dndscv.patch && \
    cd .. && \
    tar --exclude='.git' -cvzf dndscv.tar.gz dndscv

FROM debian:buster-slim AS final

ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get -y install --no-install-recommends --no-install-suggests \
        locales-all=2.28-10+deb10u4 \
        r-base=3.5.2-1 \
        r-cran-devtools=2.0.1-1 \
        r-cran-ade4=1.7-13-1 \
        r-cran-seqinr=3.4-5-2+b1 \
        r-bioc-biostrings=2.50.2-1 \
        r-bioc-genomicranges=1.34.0+dfsg-1 \
        r-bioc-rsamtools=1.34.1-1 \
        procps=2:3.3.15-2 && \
    rm -rf /var/lib/apt/lists/*

RUN --mount=type=bind,from=build,target=/build R --vanilla <<EOF
tryCatch(
{
    library(devtools)
    install.packages("/build/dndscv.tar.gz", repos=NULL, type="source")
}, error = function(e) {
    cat('Error:', e$message, '\n')
    quit(status = 1)
})
EOF

RUN --mount=type=bind,target=/files \
    mkdir -p /dndscv && \
    cp /files/dndscv.R /dndscv/

ENTRYPOINT [ "Rscript", "/dndscv/dndscv.R" ]
