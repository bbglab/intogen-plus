import os
import json
import click
import numpy as np
from functools import partial
from tqdm import tqdm
import pickle
import gzip
import csv
import itertools
from bgoncotree.main import BGOncoTree

import pandas as pd
# from pathos.multiprocessing import Pool
from multiprocessing import Pool

from utils import complementary, mut_key_generator, normalize_profile, shortkey_to_lex
from config import Values

# global variables

folder = os.environ.get("INTOGEN_DATASETS")  # './'

SITE_COUNTS_PATH = Values.SITE_COUNTS_PATH

TRI_COUNT_GENOME_PATH = Values.TRI_COUNT_GENOME_PATH
with gzip.open(TRI_COUNT_GENOME_PATH, 'rb') as f:
    tri_count_genome = json.load(f)

TRI_COUNT_EXOME_PATH = Values.TRI_COUNT_EXOME_PATH
with gzip.open(TRI_COUNT_EXOME_PATH, 'rb') as f:
    tri_count_exome = json.load(f)

SIGNATURES_PATH = Values.SIGNATURES_PATH
signatures = pd.read_csv(SIGNATURES_PATH, sep='\t', index_col=0)
signatures = {sign: signatures.loc[sign, :].values.tolist() for sign in signatures.index}


class dNdSOut:

    def __init__(self, annotmuts, genemuts):
        self.annotmuts = pd.read_csv(annotmuts, sep='\t')
        self.genemuts = pd.read_csv(genemuts, sep='\t')
        self.samples = self.annotmuts['sampleID'].unique()

    @staticmethod
    def pyr_context(row):
        """pyr-context associated to row in maf table"""

        context = row['ref3_cod'] + '>' + row['mut_cod']
        if context[1] not in list('CT'):
            context = complementary(context)
        return context

    @property
    def catalogue(self):
        """matrix with counts per sample and pyr-context"""

        df = self.annotmuts[(self.annotmuts['ref'].isin(list('ACGT'))) & (self.annotmuts['alt'].isin(list('ACGT')))]
        df['context'] = df.apply(self.pyr_context, axis=1)
        catalogue = pd.crosstab(df.sampleID, df.context)
        return catalogue

    @property
    def burden(self):
        """dict mutation count per sample"""

        dg = self.annotmuts.groupby('sampleID').count()
        burden = dict(zip(dg.index.tolist(), dg.values.flatten()))
        return burden

    def expected_syn(self, gene):
        """
        Args:
            gene: str: gene symbol
        Returns:
            expected syn mutation count at gene predicted by negative-binomial model
        """
        return self.genemuts[self.genemuts['gene_name'] == gene]['exp_syn_cv'].values[0]

    @property
    def relative_syn(self):
        """dict of proportion of syn mutations per sample"""

        df = self.annotmuts
        samples = df['sampleID'].unique()
        syn_burden = {s: len(df[(df['sampleID'] == s) & (df['impact'] == 'Synonymous')]) for s in samples}
        return {s: syn_burden[s] / sum(syn_burden.values()) for s in syn_burden}


class SiteCounts:

    def __init__(self):

        with gzip.open(SITE_COUNTS_PATH, 'rb') as f:
            self.site_counts = pickle.load(f)

    def get(self, gene, csqn_type):
        """
        dictionary: count per context
                    pyrimidine-centered lex context key -> counts
        """

        i = self.site_counts['genes'].index(gene)
        k = self.site_counts['csqn_types'].index(csqn_type)
        counts = self.site_counts['matrix'][i, :, k]
        count_dict = dict(zip(self.site_counts['contexts'], counts))
        for k in count_dict:  # e.g. count_dict['ACT>T'] -> 10
            if k[1] in list('CT'):
                count_dict[k] += count_dict[complementary(k)]
        lex_count_dict = {shortkey_to_lex(k): count_dict[k] + 1 for k in count_dict if k[1] in {'C', 'T'}}
        return lex_count_dict

    def context_per_gene(self, gene):
        """
        dictionary: count per gene-context
        """

        i = self.site_counts['genes'].index(gene)
        counts = self.site_counts['matrix'].sum(axis=2)[i, :]
        count_dict = dict(zip(self.site_counts['contexts'], counts))
        for k in count_dict:  # e.g. count_dict['ACT>T'] -> 10
            if k[1] in list('CT'):
                count_dict[k] += count_dict[complementary(k)]
        d = {shortkey_to_lex(k): count_dict[k] + 1 for k in count_dict if k[1] in {'C', 'T'}}
        d = [d[k] for k in mut_key_generator()]
        return d


def combine_signatures(weights, sig_labels, scope='exome'):
    """
    Args:
        weights: array of relative exposures
        sig_labels: list of signature labels
        scope: which region should we normalize for
    Returns:
        normalized linear combination of signatures according to weights
    """

    combined = {mk: 0 for mk in mut_key_generator()}
    for i, sign in enumerate(sig_labels):
        profile = signatures[sign]
        for j, mk in enumerate(mut_key_generator()):
            combined[mk] += weights[i] * profile[j]
    total = sum(combined.values())
    combined = {k: combined[k] / total for k in combined}
    if scope == 'exome':
        combined = normalize_profile(combined, tri_count_exome)
    elif scope == 'genome':
        combined = normalize_profile(combined, tri_count_genome)
    return combined


def genewise_run(gene, weights, annotmuts, genemuts):
    dndsout = dNdSOut(annotmuts, genemuts)
    syn_sites = SiteCounts().get(gene, 'synonymous_variant')  # counts per context
    exp = dndsout.expected_syn(gene)  # expected syn count at gene
    syn_proportion = dndsout.relative_syn  # proportion of syn mutations per sample
    expected_syn = {sample: prop * exp for sample, prop in syn_proportion.items()}
    weights_df = pd.read_csv(weights, sep='\t', index_col=0)

    # weights table after signature fitting
    # cols have signatures + SSE + mut_count

    rates = {sample: combine_signatures(weights_df.loc[sample, :].values, weights_df.columns[:-2].tolist()) for sample in dndsout.samples if sample in weights_df.index}
    mutrate = {}
    for sample in rates:
        k = sum([rates[sample][ctxt] * syn_sites[ctxt] for ctxt in mut_key_generator()])
        k = expected_syn[sample] / k
        mutrate[sample] = list(map(lambda x: k * x, [rates[sample][ctxt] for ctxt in mut_key_generator()]))
    return {gene: mutrate}


def normalization_constant(mutrate_output_folder, output_json):
    """
    computes the normalization constant for each sample
    this is a necessary first step in function normalization()

    mutrate_output_folder: folder where mutrate results are kept for each gene
    output_json: json file with dict with normalization constant per sample

    Example
    -------
    normalization_constant('./mutrate_output', './constants.json')
    """

    site_counts = SiteCounts()
    res = {}
    for fn in tqdm(os.listdir(mutrate_output_folder)):
        with open(os.path.join(mutrate_output_folder, fn), 'rt') as f:
            d = json.load(f)
        gene = next(iter(d.keys()))
        context_count = site_counts.context_per_gene(gene)
        context_count = np.array(context_count)
        for sample in d[gene]:
            arr = np.array(d[gene][sample])
            res[sample] = res.get(sample, 0) + np.dot(context_count, arr)
    with open(output_json, 'wt') as f:
        json.dump(res, f)


@click.group()
def cli():
    pass


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--mutrate_output_folder', type=click.Path(), help='results folder')
@click.option('--constants_path', type=click.Path(), help='json file for constants dict')
def normalization(mutrate_output_folder, constants_path):
    """
    computes a normalized version of mutrate, i.e.:
    conditioning on the observation of a single mutation, what is the per-site probability
    for the mutation to occur across all possible coding mutation sites

    Example
    -------
    normalization('./mutrate_output', './constants.json')
    """

    normalization_constant(mutrate_output_folder, constants_path)
    with open(constants_path, 'rt') as f:
        constants = json.load(f)
    for fn in tqdm(os.listdir(mutrate_output_folder)):
        with open(os.path.join(mutrate_output_folder, fn), 'rt') as f:
            d = json.load(f)
        gene = next(iter(d.keys()))
        norm_d = {gene: {}}
        for sample in d[gene]:
            norm_d[gene][sample] = [x / constants[sample] if constants[sample] != 0 else 0 for x in d[gene][sample]]
        with open(os.path.join(mutrate_output_folder, 'norm_' + fn), 'wt') as g:
            json.dump(norm_d, g)


@cli.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--annotmuts', type=click.Path(), help='path to dndsout$annotmuts')
@click.option('--genemuts', type=click.Path(), help='path to dndsout$genemuts')
@click.option('--weights', type=click.Path(), help='path to signature weights upon fitting')
@click.option('--drivers', type=click.Path(), help='path to the drivers file')
@click.option('--oncotree', type=click.Path(), help='path to oncotree file')
@click.option('--cohorts', type=click.Path(), help='path to cohorts file')
@click.option('--cores', default=os.cpu_count(), help='Max processes to run multiprocessing', type=click.INT)
@click.option('--output', required=True, type=click.Path(), help='file folder for outputs')
@click.option('--test', is_flag=True, help='test flag to run a reduced number of genes')
def compute_mutrate(annotmuts, genemuts, weights, drivers, oncotree, cohorts, cores, output, test):
    """
    Requirements
    ------------
    compute_mutrate.py uses the fitting obtained by deconstructSig

    Example
    -------
    python compute_mutrate.py compute-mutrate --annotmuts <annotmuts_path> --genemuts <genemuts_path> \\
                                              --weights <deconstructSigs_path> --drivers <drivers_path>
                                              -- oncotree <oncotree_path> --cohorts <stats_cohorts_path> cores <multiprocessing_cores> \\
                                              --output <output_folder>
    """
    dndsout = dNdSOut(annotmuts, genemuts)

    # Define the geneset: minimal set required to annotate MAF with expected mutations per bp

    # filter only the driver genes (using the gene symbol) from a list of cohorts and output a stingle JSON file with the genes

    # example of annotmuts value: "./annotmuts_files/CBIOP_WGS_PRAD_EURUROL_2017.annotmuts"
    file_name = annotmuts.split("/")[-1]
    cohort_str = file_name.split(".")[0]
    tree = BGOncoTree(oncotree, "\t")

    driver_genes_dict = dict()
    # drivers.tsv
    # SYMBOL	TRANSCRIPT	COHORT	CANCER_TYPE	METHODS	MUTATIONS	SAMPLES	%_SAMPLES_COHORT	QVALUE_COMBINATION	ROLE	CGC_GENE	CGC_CANCER_GENE	DOMAIN	2D_CLUSTERS	3D_CLUSTERS	EXCESS_MIS	EXCESS_NON	EXCESS_SPL
    # ABCB1	ENST00000622132	ICGC_WGS_ESAD_UK	ESAD	dndscv,cbase	16.0	14.0	0.09333333333333334	2.1848643389296567e-05	Act	False	False		87520818:87520818		0.9721893272234644	0.0	0.0
    with open(drivers, newline='') as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter="\t")
        for row in csv_reader:
            if row['COHORT'] in driver_genes_dict.keys():
                driver_genes_dict[row['COHORT']].append(row['SYMBOL'])
            else:
                array = [row['SYMBOL']]
                driver_genes_dict[row['COHORT']] = array

    tumor_type = get_tumor_type_by_cohort(cohort_str, cohorts)
    cohorts_list = get_cohorts_by_tumor_type_with_children(cohorts, tree, tumor_type)  # includes also the current cohort

    ## We want to include more driver genes. Genes that are not drivers in this cohort but in others that are related to each other via the cancer type (oncotree)
    drivers_from_cohorts = get_driver_genes_from_cohort_list(cohorts_list, driver_genes_dict)
    #print("\n list of cohorts\n", cohorts_list)
    #print("\n list of drivers from all cohorts:\n", drivers_from_cohorts)

    dg = dndsout.annotmuts[dndsout.annotmuts['mut'].isin(list('ACGT'))]
    dg = dndsout.annotmuts[dndsout.annotmuts['gene'].isin(drivers_from_cohorts)]

    gene_set = dg['gene'].unique()

    if test:
        gene_set = gene_set[:1]

    # instantiate a partial of genewise_run that encapsulates the main task
    task = partial(genewise_run, weights=weights, annotmuts=annotmuts, genemuts=genemuts)

    # prepare output dir
    if not os.path.exists(output):
        os.makedirs(output)

    # example of annotmuts value: "./annotmuts_files/CBIOP_WGS_PRAD_EURUROL_2017.annotmuts"
    output_file_name = cohort_str + ".out.json.gz"

    whole_res = dict()
    with Pool(cores) as pool:
        for res in tqdm(pool.imap(task, gene_set), total=len(gene_set)):
            whole_res.update(res)

    # dump to 1 zipped
    # json file per cohort with all genes , KEY:gene , VALUE:list of values
    with gzip.open(os.path.join(output, output_file_name), 'wt') as f_output:
        json.dump(whole_res, f_output)



def get_driver_genes_from_cohort_list(cohort_list, drivers_dict):
    list_of_drivers = []
    for cohort in cohort_list:
        if cohort in drivers_dict.keys():
            list_of_drivers.append(drivers_dict[cohort])

    flattened_list_of_drivers = list(itertools.chain(*list_of_drivers))
    return flattened_list_of_drivers


# stats_cohorts.tsv
# AGE	COHORT	MUTATIONS	PLATFORM	SAMPLES	TREATED	TYPE	CANCER_TYPE	SOURCE	LEGEND
# Adult	CBIOP_WXS_BRCA_MBCPROJECT_2018_PRY_TREAT	349	WXS	6	Treated	Primary	BRCA	CBIOP	Breast Adeno
def get_tumor_type_by_cohort(cohort, cohorts_file):
    stats_cohorts_df = pd.read_csv(cohorts_file, sep="\t")
    tumor_types = stats_cohorts_df[stats_cohorts_df["COHORT"] == cohort]["CANCER_TYPE"].unique()
    if len(tumor_types) > 0:
        return tumor_types[0]
    else:
        print("No cohort ", cohort, " found in the cohorts file")
        return None


# oncotree_boostDM.tsv
# ID	PARENT	NAMES	TAGS
# ACC	SOLID	Adrenocortical carcinoma
def get_descendant_string_list(tree, node_str):
    # given a tumor type (string) it returns a list (of strings) with all its children and children of children
    list_of_descendants = []
    dictionary_descendants = dict()
    children = tree.descendants(node_str) # it returns a list of node objects with the full tree path ex. CANCER/NON_SOLID/L/LY/DLBCL  , CANCER/NON_SOLID/L/LY/BLY  when node_str is LY (So, the parents will be removed to keep only the children)
    for child in children:
        index_of_node = str(child).find(node_str)
        break
    pos_in_string = index_of_node + len(str(node_str)) + 1
    for child in children:
        descendant = str(child)[pos_in_string:] # from string CANCER/NON_SOLID/L/LY/DLBCL I want to keep only string DLBCL (it could have more levels, from DLBCL onwards.. all children after LY)
        list_of_descendants.append(descendant)

    for strings_tree in list_of_descendants:
        nodes = strings_tree.split("/")
        for node in nodes:
            dictionary_descendants[node] = 1

    return list(dictionary_descendants.keys())


def get_cohorts_by_tumor_type_with_children(cohorts_file, tree, ttype):
    # given a tumor type (string) it returns all cohorts found with the tumor type in the stats_cohort.tsv PLUS the cohorts of the children of the tumor
    stats_cohorts = pd.read_csv("stats_cohorts.tsv", sep="\t")
    children_ttypes = get_descendant_string_list(tree, ttype)
    cohorts = stats_cohorts[stats_cohorts["CANCER_TYPE"] == ttype]["COHORT"].unique()
    all_cohorts = cohorts
    if children_ttypes is None:
        return all_cohorts
    else:  # if there are children
        for child_tumor in children_ttypes:
            child_cohorts = get_cohorts_by_tumor_type(cohorts_file , child_tumor)
            all_cohorts = np.concatenate((child_cohorts, all_cohorts))

        return all_cohorts


def get_cohorts_by_tumor_type(cohorts_file, ttype):
    # given a tumor type it returns all cohorts found with the tumor type in the stats_cohort.tsv PLUS the cohorts of the children of the tumor
    stats_cohorts = pd.read_csv(cohorts_file, sep="\t")
    cohorts = stats_cohorts[stats_cohorts["CANCER_TYPE"] == ttype]["COHORT"].unique()
    return cohorts



if __name__ == '__main__':
    cli()

    # tests
    # -----
    # normalization_constant('./mutrate_output', './constants.json')
    # normalization('./mutrate_output', './normalizing_constant.json')

# HOW TO DO A TEST RUN:
# export INTOGEN_DATASETS="./hg38_vep92_develop"
# python compute_mutrate.py compute-mutrate --annotmuts ./annotmuts_files/CBIOP_WGS_PRAD_EURUROL_2017.annotmuts --genemuts  ./genemuts/CBIOP_WGS_PRAD_EURUROL_2017.genemuts --weights ./CBIOP_WGS_PRAD_EURUROL_2017.out  --cores 3 --output ./CBIOP_WGS_PRAD_EURUROL_2017_outs --drivers drivers.tsv --cohorts stats_cohorts.tsv --oncotree oncotree_boostDM.tsv