FROM python:3.9

RUN mkdir -p /combination
COPY LICENSE.txt /combination/
COPY MANIFEST.in /combination/
COPY README.rst /combination/
COPY requirements.txt /combination/
COPY setup.py /combination/
COPY intogen_combination /combination/intogen_combination

RUN pip install --no-cache-dir /combination

ENTRYPOINT [ "/usr/local/bin/intogen-combine" ]
