import gzip
import dill as pickle
import pandas as pd

from os import path


class Reader(object):

    """
    Reads paths from configuration file and keeps them as attributes
    It has methods to instantiate inputs from those paths
    """

    def __init__(self):
        config = {
            'base': {
                'site_counts_file': path.join(os.environ['INTOGEN_DATASETS'], 'mutrate', 'site_counts.pickle.gz'),
                'triplets': path.join(os.environ['INTOGEN_DATASETS'], 'mutrate','tnt_counts_dict.pickle.gz'),
                'hugo_ensembl': path.join(os.environ['INTOGEN_DATASETS'], 'mutrate','HGNC_ENSEMBL_dict.pickle'),
                'ensembl_hugo': path.join(os.environ['INTOGEN_DATASETS'], 'mutrate','ENSEMBL_HGNC_dict_70.pickle')
            },
            'sigfit': {
                'cosmic': path.join(os.environ['INTOGEN_DATASETS'], 'mutrate','cosmic_signatures.pickle.gz'),
                'cosmic_exome': path.join(os.environ['INTOGEN_DATASETS'], 'mutrate','cosmic_exome.tsv')
            }
        }

        self.site_counts_path = config['base']['site_counts_file']
        self.triplets_path = config['base']['triplets']
        self.hugo_ensembl_path = config['base']['hugo_ensembl']
        self.ensembl_hugo_path = config['base']['ensembl_hugo']
        self.cosmic_path = config['sigfit']['cosmic']
        self.cosmic_exome_path = config['sigfit']['cosmic_exome']

    def read_site_counts(self):
        with gzip.open(self.site_counts_path, 'rb') as f_input:
            site_counts = pickle.load(f_input)
        return site_counts

    def read_triplets(self, key):
        with gzip.open(self.triplets_path, 'rb') as f_input:
            triplets = pickle.load(f_input)
        return triplets[key]

    def read_cosmic(self):
        with gzip.open(self.cosmic_path, 'rb') as f_input:
            cosmic_sigs = pickle.load(f_input)
        return cosmic_sigs

    @staticmethod
    def read_exposures(mut_catalogue_path, sigfit_exposures_path):
        df_mean = pd.read_csv(sigfit_exposures_path, sep='\t')
        dg_mutations = pd.read_csv(mut_catalogue_path, sep='\t', index_col=0)
        dg_mutations['total'] = dg_mutations.sum(axis=1)
        df_mean.index = dg_mutations.index
        explained_muts_mean = df_mean.values.T * dg_mutations['total'].values
        dh_mean = pd.DataFrame(explained_muts_mean.T, index=df_mean.index, columns=df_mean.columns)
        return dh_mean

    def read_gene_dicts(self):
        with open(self.hugo_ensembl_path, 'rb') as f_input:
            he_dict = pickle.load(f_input)
        with open(self.ensembl_hugo_path, 'rb') as f_input:
            eh_dict = pickle.load(f_input)
        return he_dict, eh_dict
