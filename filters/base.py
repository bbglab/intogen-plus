from bgparsers import readers


class Filter:

    def run(self, group_key, group_data):
        """
        :param group_key: Group key identifications
        :param group_data: Group data parameters
        :return: A rows generator as a dictionary
        """

        # This method must be implemented at subclass level
        raise NotImplementedError()


class VariantsReader(Filter):

    def run(self, group_key, group_data):
        return readers.variants(group_data)
