FROM debian:buster-slim

ENV LC_ALL=C.UTF-8

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        python3-pandas python3-click r-bioc-bsgenome r-cran-gplots procps && \
    rm -rf /var/lib/apt/lists/*

SHELL ["/bin/bash", "-o", "pipefail", "-c"]
RUN echo "install.packages(\"BiocManager\", repos=\"https://cran.r-project.org\")" | R --no-save
RUN echo "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")" | R --no-save
RUN echo "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")" | R --no-save
RUN echo "install.packages(\"deconstructSigs\", repos=\"https://cran.r-project.org\")" | R --no-save

RUN mkdir -p /deconstructsig
COPY run_deconstruct.py /deconstructsig/
COPY deconstructSigs.r /deconstructsig/
COPY signature_assignment.py /deconstructsig/
COPY output_pass_drivers_01.csv /deconstructsig/
RUN chmod -R a+rx /deconstructsig/

ENTRYPOINT [ "/usr/bin/python3", "/deconstructsig/run_deconstruct.py" ]
