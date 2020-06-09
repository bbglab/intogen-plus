import csv
import gzip


def _open(file):
    if file.endswith('.gz'):
        return gzip.open(file, 'rt')
    else:
        return open(file, 'r')


class TSVReader:
    def __init__(self, file):
        self._filename = file
        self._f = None

    def __enter__(self):
        self._f = _open(self._filename)
        reader = csv.DictReader(self._f, delimiter='\t')
        return reader

    def __exit__(self, exc_type, exc_value, traceback):
        self._f.close()

    def __iter__(self):
        with self as reader:
            for row in reader:
                yield row
