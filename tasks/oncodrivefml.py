import os
import csv
import gzip
import shutil

from os import path
from .base import Task
from oncodrivefml.main import OncodriveFML
from oncodrivefml.config import load_configuration


class OncodriveFmlTask(Task):

    KEY = 'oncodrivefml'

    def __init__(self, output_folder, config):
        super().__init__(output_folder, config)

        self.name = None
        self.config = load_configuration(config[OncodriveFmlTask.KEY]['config_file'], override={'settings': {'cores': 2}})
        self.elements_file = self.config['elements']

        self.in_file = None
        self.in_skip = False
        self.in_fd = None
        self.in_writer = None

        self.out_file = None
        self.output_folder = path.join(output_folder, "oncodrivefml")
        os.makedirs(self.output_folder, exist_ok=True)

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

            # Check if whole genome
            if "_WGS_" in os.path.basename(self.in_file):
                self.config['signature']['normalize_by_sites'] = 'whole_genome'

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

    def clean(self):

        if path.exists(self.out_file):
            os.unlink(self.out_file)

        if path.exists(self.in_file):
            os.unlink(self.in_file)

        tmp_folder = self.out_file + ".tmp"
        if path.exists(tmp_folder):
            shutil.rmtree(tmp_folder)

