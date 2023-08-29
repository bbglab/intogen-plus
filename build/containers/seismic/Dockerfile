FROM rocker/r-ver:4.1

ENV LC_ALL=C.UTF-8
ENV R_BASE_VERSION=4.1.0

RUN mkdir /seismic/

COPY install_R_dependencies.R /seismic/
COPY SEISMIC.R /seismic/
COPY data/ /seismic/data/
COPY config/ /seismic/config/

RUN chmod -R a+rx /seismic/

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        r-base-html \
        r-doc-html \
        r-cran-devtools \
        r-bioc-biostrings \
        r-bioc-genomicranges \
        r-bioc-rsamtools \
        r-cran-ade4 \
        r-cran-seqinr \
        procps apt-utils \
        locales-all \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev  \
        libfribidi-dev  \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libcairo2-dev \
        libxt-dev \
        libmagick++-dev

RUN Rscript /seismic/install_R_dependencies.R


# ENTRYPOINT ["Rscript", "/seismic/SEISMIC.R"]

# CMD "$@"