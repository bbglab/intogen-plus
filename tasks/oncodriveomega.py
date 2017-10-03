import csv
import gzip
import os

from os import path
from oncodriveomega.config import load_configuration
from oncodriveomega.main import OncodriveOmega


class OncodriveOmegaTask:

    def __init__(self, output_folder, config_file):

        self.name = None
        self.config_file = config_file
        # TODO add expression filter in config file
        self.expression_filter = None

        self.in_fd = None
        self.in_writer = None
        self.in_file = None
        self.in_skip = False

        self.out_file = None
        self.output_folder = path.join(output_folder, 'oncodriveomega')
        os.makedirs(self.output_folder, exist_ok=True)

    def __repr__(self):
        return "OncodriveOmega '{}'".format(self.name)

    def init(self, name):
        self.name = name
        self.in_file = path.join(self.output_folder, "{}.in.gz".format(name))
        self.in_skip = path.exists(self.in_file)
        self.out_file = path.join(self.output_folder, "{}.out.gz".format(name))

    def input_start(self):

        if not self.in_skip:
            self.in_fd = gzip.open(self.in_file, 'wt')
            self.in_writer = csv.writer(self.in_fd, delimiter='\t')
            self.in_writer.writerow(['SAMPLE', 'CHR', 'POS', 'REF', 'ALT', 'GENE', 'CSQN'])

    def input_write(self, _, value):

        if not self.in_skip:

            # SAMPLE CHR POS REF ALT GENE CSQN
            if value['CANONICAL'] == 'YES':
                identifier, sample, ref, alt = value['#Uploaded_variation'].split('__')
                chromosome, position = value['Location'].split(':')
                consequence = value['Consequence'].split(',')[0]

                self.in_writer.writerow([
                    sample,
                    chromosome,
                    position,
                    ref,
                    alt,
                    value['Gene'],
                    consequence
                ])

    def input_end(self):

        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Load config
        config = load_configuration(self.config_file)

        signature_file = path.join(config['datasets']['signatures_folder'],
                                   "{}.txt_signature_full.pickle.gz".format(self.name))

        if not path.exists(signature_file):
            signature_file = None

        covariates_file = path.join(config['datasets']['covariates_folder'], "{}.covariates_ensembl.tsv".format(self.name))
        if not path.exists(covariates_file):
            covariates_file = path.join(config['datasets']['covariates_folder'], "GENERIC.covariates_ensembl.tsv")

        analysis = OncodriveOmega(self.in_file, signature_file, covariates_file, config=config, cores=4)
        analysis.prepare_mutations(cache=False)
        analysis.prepare_regression(cache=False)

        # Load all mutated elements
        with gzip.open(self.in_file, 'rt') as fd:
            reader = csv.reader(fd, delimiter='\t')
            next(reader)
            elements = set([e[5].strip() for e in reader])

        # Genes in regression
        genes_in_regression = set(analysis.regression_dataframe.index)
        elements = elements.intersection(genes_in_regression)

        # Genes with expression filter
        if self.expression_filter is not None:
            with open(self.expression_filter, 'rt') as fd:
                expressed = set([e[0].strip() for e in csv.reader(fd, delimiter='\t')])
            elements = elements.intersection(expressed)

        results = analysis.run(elements)

        analysis.wrap_results(results, self.out_file)

    def clean(self):

        if path.exists(self.out_file):
            os.unlink(self.out_file)
