FROM zhouwanding/transvar:2.5.9

ENV LC_ALL=C.UTF-8
ENV TRANSVAR_CFG=/data/transvar.cfg
ENV TRANSVAR_DOWNLOAD_DIR=/data

RUN rm -rf /anno && \
    mkdir -p /data

ENTRYPOINT [ "/usr/local/bin/transvar" ]
