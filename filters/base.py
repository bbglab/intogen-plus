import csv
import gzip

from bgparsers import readers


class Filter:

    def __init__(self, parent):
        self.parent = parent
        self.stats = {}

    def run(self, group_key, group_data):
        """
        :param group_key: Group key identifications
        :param group_data: Group data parameters
        :return: A rows generator as a dictionary
        """

        # This method must be implemented at subclass level
        raise NotImplementedError()


class VariantsReader(Filter):

    def __init__(self):
        super().__init__(None)

    def run(self, group_key, group_data):
        return readers.variants(group_data)


class TSVReader(Filter):

    def __init__(self):
        super().__init__(None)

    def run(self, group_key, tsv_file):
        with gzip.open(tsv_file, "rt") as fd:
            reader = csv.DictReader(fd, delimiter='\t')
            for row in reader:
                yield row
