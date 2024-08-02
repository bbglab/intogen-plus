FROM debian:buster-slim

ENV LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get -y install --no-install-recommends \
        locales-all r-base r-cran-devtools r-bioc-biostrings r-bioc-genomicranges \
        r-bioc-rsamtools r-cran-ade4 r-cran-seqinr procps && \
    rm -rf /var/lib/apt/lists/*

RUN mkdir -p /dndscv
COPY dndscv.R /dndscv/
# COPY dndscv-ae10c76.tar.gz /dndscv/
RUN chmod -R a+rx /dndscv/

# # RUN echo "library(devtools); install_github('im3sanger/dndscv')" | R --no-save
# RUN echo "library(devtools)" | R --no-save
# RUN R CMD INSTALL "/dndscv/dndscv-ae10c76.tar.gz"

ENTRYPOINT [ "Rscript", "/dndscv/dndscv.R" ]
