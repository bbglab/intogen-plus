FROM python:3

RUN pip install --no-cache-dir open-variant
# RUN pip install --no-cache-dir https://bitbucket.org/bgframework/bgparsers/get/master.tar.gz

RUN mkdir -p /core
COPY LICENSE.txt /core/
COPY MANIFEST.in /core/
COPY README.rst /core/
COPY requirements.txt /core/
COPY setup.py /core/
COPY intogen_core /core/intogen_core

RUN pip install --no-cache-dir /core

COPY get_field.sh /usr/local/bin/
RUN chmod a+rx /usr/local/bin/get_field.sh
