import os
import csv
import gzip
import subprocess

from os import path
from .base import Task


class EDriverTask(Task):

    KEY = 'edriver'

    def __init__(self, output_folder, config):

        super().__init__(output_folder, config)

        self.name = None
        self.in_fd = None
        self.in_writer = None
        self.in_file = None
        self.in_skip = False

        self.out_file = None
        self.output_folder = path.join(output_folder, self.KEY)
        os.makedirs(self.output_folder, exist_ok=True)

    def input_start(self):

        if not self.in_skip:
            self.in_fd = gzip.open(self.in_file, 'wt')
            self.in_writer = csv.writer(self.in_fd, delimiter='\t', lineterminator='\n')

    def input_write(self, _, value):

        if not self.in_skip:

            ensp, position = value['ENSP'], value['Protein_position']
            if ensp != "-" and position != "-":
                identifier, sample, ref, alt = value['#Uploaded_variation'].split('__')

                # ENSP Protein_position sample tissue
                self.in_writer.writerow([
                    ensp,
                    position,
                    sample,
                    "unknown"
                ])

    def input_end(self):
        if not self.in_skip:
            self.in_fd.close()
            self.in_writer = None
            self.in_fd = None

    def run(self):

        # Run vep
        if not path.exists(self.out_file):
            cmd = "mkdir -p {2}/tmp_{3} && " \
                  "zcat {1} > {2}/tmp_{3}/input.txt && " \
                  "singularity run {0}/edriver.simg 1 {2}/tmp_{3}/input.txt {0}/iur_and_pfam_ensembl_v85.txt {0}/ensembl_v85_seqs_parsed_header.txt {2}/tmp_{3}/ed &&" \
                  "bash -c 'cat <(head -n1 {2}/tmp_{3}/ed_unknown_with_corrected_p_values.txt) <(grep unknown {2}/tmp_{3}/ed_unknown_with_corrected_p_values.txt)' | cut -f2- | gzip > {2}/{3}.out.gz &&" \
                  "rm -rf {2}/tmp_{3}".format(
                os.environ['EDRIVER_DATA'], self.in_file, os.path.abspath(self.output_folder), self.name
            )

            stdout = subprocess.check_output(cmd, shell=True)
            print(stdout.decode())

            with gzip.open("{}/mart_export_85.txt.gz".format(os.environ['EDRIVER_DATA']), 'rt') as fd:
                prot_to_gene = {v['Ensembl Protein ID']: v['Ensembl Gene ID'] for v in csv.DictReader(fd, delimiter='\t')}

            with gzip.open(os.path.abspath(self.out_file), 'rt') as fd:
                results = {}
                for row in csv.DictReader(fd, delimiter='\t'):
                    if row['Protein'] in prot_to_gene:
                        gene = prot_to_gene[row['Protein']]
                        v = results.get(gene, (2.0, 2.0))
                        if float(row['p']) < v[0]:
                            results[gene] = (float(row['p']), float(row['q']))

            with gzip.open("{}/{}.genes.out.gz".format(self.output_folder, self.name), 'wt') as fd:
                writer = csv.writer(fd, delimiter='\t')
                writer.writerow(['GENE', 'PVALUE', 'QVALUE'])
                for gene, vals in sorted(results.items(), key=lambda v: v[1][0]):
                    writer.writerow([gene, vals[0], vals[1]])

        return self.out_file
