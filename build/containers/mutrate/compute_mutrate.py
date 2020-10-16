import os
import json
import click
from functools import partial
from tqdm import tqdm
import pickle
import gzip
import csv
import itertools

from collections import defaultdict

import pandas as pd
from multiprocessing import Pool
from bgoncotree.main import BGOncoTree

from utils import complementary, mut_key_generator, normalize_profile, shortkey_to_lex
from config import Values

# global variables

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

    rates = {sample: combine_signatures(weights_df.loc[sample, :].values, weights_df.columns[:-2].tolist()) for sample
             in dndsout.samples if sample in weights_df.index}
    mutrate = {}
    for sample in rates:
        k = sum([rates[sample][ctxt] * syn_sites[ctxt] for ctxt in mut_key_generator()])
        k = expected_syn[sample] / k
        mutrate[sample] = list(map(lambda x: k * x, [rates[sample][ctxt] for ctxt in mut_key_generator()]))
    return {gene: mutrate}


@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--annotmuts', type=click.Path(), help='path to dndsout$annotmuts')
@click.option('--genemuts', type=click.Path(), help='path to dndsout$genemuts')
@click.option('--weights', type=click.Path(), help='path to signature weights upon fitting')
@click.option('--drivers', type=click.Path(), help='path to the drivers file')
@click.option('--oncotree', type=click.Path(), help='path to oncotree file')
@click.option('--cohorts', type=click.Path(), help='path to cohorts file')
@click.option('--cores', default=os.cpu_count(), help='Max processes to run multiprocessing', type=click.INT)
@click.option('--test', is_flag=True, help='test flag to run a reduced number of genes')
def cli(annotmuts, genemuts, weights, drivers, oncotree, cohorts, cores, test):
    """
    Requirements
    ------------
    compute_mutrate.py uses the fitting obtained by deconstructSig

    Example
    -------
    python compute_mutrate.py  --annotmuts <annotmuts_path> --genemuts <genemuts_path> \\
                                --weights <deconstructSigs_path> --drivers <drivers_path>
                                -- oncotree <oncotree_path> --cohorts <stats_cohorts_path> cores <multiprocessing_cores> \\
    """
    dndsout = dNdSOut(annotmuts, genemuts)

    # Define the geneset: minimal set required to annotate MAF with expected mutations per bp

    # filter only the driver genes (using the gene symbol) from a list of cohorts and output a stingle JSON file with the genes

    # example of annotmuts value: "./annotmuts_files/CBIOP_WGS_PRAD_EURUROL_2017.annotmuts"
    file_name = annotmuts.split("/")[-1]
    cohort_str = file_name.split(".")[0]
    tree = BGOncoTree(oncotree, "\t")

    # drivers.tsv
    # SYMBOL	TRANSCRIPT	COHORT	CANCER_TYPE	METHODS	MUTATIONS	SAMPLES	%_SAMPLES_COHORT	QVALUE_COMBINATION	ROLE	CGC_GENE	CGC_CANCER_GENE	DOMAIN	2D_CLUSTERS	3D_CLUSTERS	EXCESS_MIS	EXCESS_NON	EXCESS_SPL
    # ABCB1	ENST00000622132	ICGC_WGS_ESAD_UK	ESAD	dndscv,cbase	16.0	14.0	0.09333333333333334	2.1848643389296567e-05	Act	False	False		87520818:87520818		0.9721893272234644	0.0	0.0
    driver_genes_dict = defaultdict(list)
    with open("drivers.tsv", newline='') as csv_file:
        csv_reader = csv.DictReader(csv_file, delimiter="\t")
        for row in csv_reader:
            driver_genes_dict[row['COHORT']].append(row['SYMBOL'])

    cohorts_df = pd.read_csv(cohorts, sep="\t")
    tumor_type = get_tumor_type_by_cohort(cohort_str, cohorts_df)
    cohorts_list = get_cohorts_by_tumor_type_with_children(cohorts_df, tree,tumor_type)  # includes also the current cohort

    # We want to include more driver genes. Genes that are not drivers in this cohort but in others that are related to each other via the cancer type (oncotree)
    drivers_from_cohorts = get_driver_genes_from_cohort_list(cohorts_list, driver_genes_dict)

    dg = dndsout.annotmuts[dndsout.annotmuts['mut'].isin(list('ACGT'))]
    dg = dndsout.annotmuts[dndsout.annotmuts['gene'].isin(drivers_from_cohorts)]

    gene_set = dg['gene'].unique()

    if test:
        gene_set = gene_set[:1]

    # instantiate a partial of genewise_run that encapsulates the main task
    task = partial(genewise_run, weights=weights, annotmuts=annotmuts, genemuts=genemuts)

    # example of annotmuts value: "./annotmuts_files/CBIOP_WGS_PRAD_EURUROL_2017.annotmuts"
    output_file_name = cohort_str + "_mutrate.out.json.gz"

    whole_res = dict()
    with Pool(cores) as pool:
        for res in tqdm(pool.imap(task, gene_set), total=len(gene_set)):
            whole_res.update(res)

    # dump to 1 zipped
    # json file per cohort with all genes , KEY:gene , VALUE:list of values
    with gzip.open(output_file_name, 'wt') as f_output:
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
def get_tumor_type_by_cohort(cohort, cohorts_df):
    tumor_types = cohorts_df[cohorts_df["COHORT"] == cohort]["CANCER_TYPE"].unique()
    if len(tumor_types) > 0:
        if len(tumor_types) > 1:
            raise NameError("in the cohorts file found more than 1 tumor type for cohort ", cohort, tumor_types)
        return tumor_types[0]
    else:
        raise NameError("No cohort ", cohort, "found in the cohorts file")


# oncotree_boostDM.tsv
# ID	PARENT	NAMES	TAGS
# ACC	SOLID	Adrenocortical carcinoma
def get_descendant_ttype_list(tree, ttype):
    """
    given a tumor type (string) it returns a list (of strings) with all its children and children of children
    """
    list_of_descendants = [ttype]
    children = tree.descendants(ttype)  # it returns a list of node objects, doesn't include the parent ttype
    for child in children:
        if child.id != ttype:
            list_of_descendants.append(child.id)

    return list_of_descendants


def get_cohorts_by_tumor_type_with_children(cohorts_df, tree, ttype):
    """
    given a tumor type (string) it returns all cohorts found with the tumor type in the stats_cohort.tsv PLUS the cohorts of the children of the tumor
    """
    children_ttypes = get_descendant_ttype_list(tree, ttype)
    all_cohorts = cohorts_df[cohorts_df["CANCER_TYPE"].isin(children_ttypes)]["COHORT"].unique()
    return all_cohorts

if __name__ == '__main__':
    cli()

# HOW TO DO A TEST RUN:
# export INTOGEN_DATASETS="./hg38_vep92_develop"
# python compute_mutrate.py --annotmuts ./annotmuts_files/CBIOP_WGS_PRAD_EURUROL_2017.annotmuts --genemuts  ./genemuts/CBIOP_WGS_PRAD_EURUROL_2017.genemuts --weights ./CBIOP_WGS_PRAD_EURUROL_2017.out  --cores 3  --drivers drivers.tsv --cohorts stats_cohorts.tsv --oncotree oncotree_boostDM.tsv
