FROM debian:buster-slim

ENV LC_ALL=C.UTF-8
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends --no-install-suggests \
        python3.7-minimal=3.7.3-2+deb10u7 \
        python3-pandas=0.23.3+dfsg-3 \
        python3-click=7.0-1 \
        r-bioc-bsgenome=1.50.0-1 \
        r-cran-gplots=3.0.1.1-1 \
        procps=2:3.3.15-2 && \
    rm -rf /var/lib/apt/lists/*

RUN R --vanilla <<EOF
tryCatch(
{
    install.packages("BiocManager", repos="https://cran.r-project.org")
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
    install.packages("https://cran.r-project.org/src/contrib/Archive/deconstructSigs/deconstructSigs_1.8.0.tar.gz", repos=NULL, type="source")
}, error = function(e) {
    cat('Error:', e$message, '\n')
    quit(status = 1)
})
EOF

RUN --mount=type=bind,target=/files \
    mkdir -p /deconstructsig && \
    cd /files && \
    cp run_deconstruct.py deconstructSigs.r signature_assignment.py /deconstructsig/ && \
    chmod -R a+rx /deconstructsig/ && \
    cp output_pass_drivers_01.csv /deconstructsig/

ENTRYPOINT [ "/usr/bin/python3", "/deconstructsig/run_deconstruct.py" ]
