FROM python:3.6

ENV BBGLAB_HOME=/oncodrivefml

RUN mkdir /oncodrivefml
COPY oncodrivefml_v2.conf /oncodrivefml
RUN chmod -R a+rx /oncodrivefml

RUN pip --no-cache-dir install oncodrivefml

WORKDIR /oncodrivefml

ENTRYPOINT [ "/usr/local/bin/oncodrivefml" ]
