from os import path


class Task:

    def __init__(self, output_folder, config):

        self.name = None
        self.in_file = None
        self.process = None
        self.in_skip = False
        self.out_file = None

        self.output_folder = output_folder
        self.config = config

    def init(self, name):

        if path.isabs(name):
            self.in_file = name
            self.name = path.basename(name).split('.')[0]
        else:
            self.name = name
            self.in_file = path.join(self.output_folder, "{}.in.gz".format(name))

        self.process = None
        self.in_skip = path.exists(self.in_file)
        self.out_file = path.join(self.output_folder, "{}.out.gz".format(self.name))

    def __repr__(self):
        return "{} '{}'".format(self.KEY, self.name)
