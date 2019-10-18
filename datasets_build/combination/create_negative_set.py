"""
Creates a negative set of genes. A gene is not considered cancer-related if:
    * is a very long gene. List taken from Przytycki and Singh, 2017
    * is an olfactory receptor. List taken from the HORDE database
    * is not expressed. Data taken from PanCanAtlas open-version

Regarding the expression of the genes:
    * samples are considered not expressed if their log2(RSEM) <= 0
    * genes are considered not expressed if 80% of the samples are not expressed
"""

# Import modules
import os
import random
from glob import glob

import click
import pandas as pd
import bgdata


class NegativeSet:

    def __init__(self, olfactory_receptors):  # path_expression_data):
        """Initialize the class
        :param olfactory_receptors: path, file with a list of olfactory receptors
        :return: None
        """
        df_olfactory_receptors = pd.read_csv(olfactory_receptors, sep='\t', header=0)
        self.olfactory_receptors = set(df_olfactory_receptors['Symbol'].tolist())
        self.path_expression_data = bgdata.get('intogen/expression/tcga_pancanatlas')  # path_expression_data
        self.negative_set, self.not_expressed = self.create_negative_set()

    def create_negative_set(self):
        """Create a negative set of genes
        :return: dictionary
        """
        results = {}
        long_genes = set([
            'HMCN1', 'TTN', 'OBSCN', 'GPR98',  'RYR2', 'RYR3'
        ])
        not_expressed_pancancer = {}
        # for input_file in glob(os.path.join(self.path_expression_data, '*.txt')):
        #     if 'samplecorr' in input_file:
        #         continue
        data = pd.read_csv(self.path_expression_data, header=0, sep='\t', compression='infer')
        tumors = data.groupby(by='TUMOR_TYPE')
        for tumor, df in tumors:
            # tumor = os.path.splitext(os.path.basename(input_file))[0]
            # df = pd.read_csv(input_file, header=0, sep='\t', compression='infer')
            try:
                values = df.groupby(by='GENE').apply(lambda x: x['log2(RSEM)'].tolist()).reset_index()
            except KeyError as e:
                print('Skipping {}'.format(tumor))
                continue
            values.rename(columns={0: 'log2(RSEM)'}, inplace=True)
            values['NOT_EXPRESSED_RATIO'] = values['log2(RSEM)'].map(
                lambda x: len([i for i in x if float(i) <= 0]) / len(x)
            )
            not_expressed = values[values['NOT_EXPRESSED_RATIO'] >= 0.8].copy()
            removed = set(not_expressed['GENE'].tolist()) | long_genes | self.olfactory_receptors
            print('{}: removed {} (total: {}, not expressed: {}'.format(
                tumor, len(removed), len(values), len(not_expressed)
            ))
            not_expressed_pancancer[tumor] = not_expressed
            results[tumor] = removed
        # Pancancer
        removed = long_genes | self.olfactory_receptors
        print('{}: removed {} (total: {}, not expressed: {}'.format(
            'PANCANCER', len(removed), len(values), 0
        ))
        results['PANCANCER'] = removed

        return results, not_expressed_pancancer

    def save(self, output, output_expression):
        """Save the negative sets to an output file
        :param output: path, path of the file to create
        :param output_expression: path, path of the output file with the list of not-expressed genes
        :return: None
        """
        with open(output, 'w') as out:
            for tumor, removed in self.negative_set.items():
                out.write('{}\t{}\n'.format(tumor, ','.join(removed)))

        # Write a list of not expressed genes
        if output_expression is not None:
            with open(output_expression, 'w') as out:
                for tumor, df in self.not_expressed.items():
                    out.write('{}\t{}\n'.format(tumor, ','.join(df['GENE'].tolist())))
                df = pd.concat(self.not_expressed.values())
                gene_list = list(set(df['GENE'].tolist()))
                out.write('{}\t{}\n'.format('PANCANCER', ','.join(gene_list)))

    @staticmethod
    def shuffle(x):
        """Add variance, between -0.05 and 0.05, to a number
        :param x: int, number
        :return: float
        """
        return x + random.uniform(-0.05, 0.05)

    def read_file(self, input_file):
        """Read a precalculad negative set file
        :param input_file: path, path of the file with the precalculated negative sets
        :return: dictionary
        """
        results = {}
        with open(input_file, 'r') as fd:
            for line in fd:
                tumor, genes = line.strip().split('\t')
                results[tumor] = genes.split(',')
        return results


@click.command()
@click.option('--olfactory_receptors', 'olfactory_receptors', help='Path to the olfactory receptors')
# @click.option('--path_expression_data', 'path_expression_data', help='Path to the expression data')
@click.option('--output_total', 'output_total', help='Output file of the total file')
@click.option('--output_non_expressed', 'output_non_expressed', help='Output file of the non-expressed genes')
def cmdline(olfactory_receptors, output_total, output_non_expressed):  # path_expression_data
    negative_set = NegativeSet(olfactory_receptors=olfactory_receptors)  # path_expression_data=path_expression_data)
    negative_set.save(output=output_total, output_expression=output_non_expressed)


if __name__ == "__main__":
    cmdline()
