import os
import csv
import gzip


from os import path

import shutil
from oncodrivefml.main import OncodriveFML
from oncodrivefml.config import load_configuration


class OncodriveFmlTask:

    ID = "oncodrivefml"

    def __init__(self, name, output_folder, config_file):

        self.name = name
        self.config = load_configuration(config_file, override={'settings': {'cores': os.cpu_count()}})
        self.elements_file = self.config['elements']

        os.makedirs(path.join(output_folder, OncodriveFmlTask.ID), exist_ok=True)
        self.in_file = path.join(output_folder, OncodriveFmlTask.ID, "{}.in.gz".format(name))
        self.in_skip = path.exists(self.in_file)
        self.in_fd = None
        self.in_writer = None

        self.out_file = path.join(output_folder, OncodriveFmlTask.ID, "{}.out.gz".format(name))

    def __repr__(self):
        return "OncodriveFML '{}'".format(self.name)

    def input_start(self):
        if not self.in_skip:
            self.in_fd = gzip.open(self.in_file, 'wt')
            self.in_writer = csv.writer(self.in_fd, delimiter='\t')
            self.in_writer.writerow(['CHROMOSOME', 'POSITION', 'REF', 'ALT', 'SAMPLE', 'ID'])

    def input_write(self, identifier, mut):
        if not self.in_skip:
            self.in_writer.writerow([
                mut['CHROMOSOME'],
                mut['POSITION'],
                mut['REF'],
                mut['ALT'],
                mut['SAMPLE'],
                identifier
            ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run OncodriveFML
        if not path.exists(self.out_file):

            tmp_folder = self.out_file + ".tmp"
            os.makedirs(tmp_folder, exist_ok=True)
            analysis = OncodriveFML(mutations_file=self.in_file,
                                    elements_file=self.elements_file,
                                    output_folder=tmp_folder,
                                    config=self.config,
                                    blacklist=None,
                                    generate_pickle=False)

            analysis.run()

            with open(analysis.output_file_prefix + ".tsv", "rb") as f_in, gzip.open(self.out_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            shutil.rmtree(tmp_folder)

        return self.out_file
