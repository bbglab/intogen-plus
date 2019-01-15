import utils
import pandas as pd
import os
from tqdm import tqdm
from pathos.multiprocessing import Pool
from functools import partial
from reader import Reader
import click
import json


def cosmic_sigfit_format(cosmic_exome_scope):
    index = ['Signature {}'.format(i+1) for i in range(30)]
    columns = [utils.lex_to_sigfit(mk) for mk in utils.mut_key_generator()]
    df = pd.DataFrame({}, columns=columns, index=index)
    for mk in utils.mut_key_generator():
        for ind in df.index:
            context = utils.lex_to_sigfit(mk)
            df.loc[ind, context] = cosmic_exome_scope[ind][mk]
    return df


def retrieve_counts(site_counts, gene_id, csqn_type='synonymous_variant'):
    """
    Args:
        site_counts: instance with counts per site by context and gene
        gene_id: str: gene ID
        csqn_type: str: consequence type
    Return:

    """
    i = site_counts['gene_index'].index(gene_id)
    k = site_counts['csqn_index'].index(csqn_type)
    counts = site_counts['matrix'][i, :, k]
    count_dict = dict(zip(site_counts['key_index'], counts))
    for k in count_dict:
        if k[1] in {'C', 'T'}:
            count_dict[k] += count_dict[utils.complementary(k)]
    return {utils.shortkey_to_lex(k): count_dict[k] + 1 for k in count_dict if k[1] in {'C', 'T'}}


def compute_rel_syn_burden(df):
    """
    Args:
        df:  dndsout$annotmuts format
    Return:
        dict: sample -> proportion of cohort syn mutations in sample
    """
    syn_burden = {sample: len(df[(df['sampleID'] == sample) & (df['impact'] == 'Synonymous')])\
                  for sample in df['sampleID'].unique()}
    total = sum(syn_burden.values())
    return {sample: syn_burden[sample] / total for sample in syn_burden}


def create_mut_catalogue(df):
    """
    Args:
        df: dataframe: dndsout$annotmuts format
    Return:
        mutational catalogue in sigfit format
    """
    mut_catalogue = utils.maf_to_makeup(df)
    index = [sample for sample in mut_catalogue]
    columns = [utils.lex_to_sigfit(mk) for mk in utils.mut_key_generator()]
    dg = pd.DataFrame({}, columns=columns, index=index)
    for mk in utils.mut_key_generator():
        for ind in dg.index:
            context = utils.lex_to_sigfit(mk)
            dg.loc[ind, context] = mut_catalogue[ind][mk]
    return dg


def switch_scope(mut_catalogue, curr_triplets, new_triplets):
    """
    Args:
        mut_catalogue: mutational catalogue
        curr_triplets: triplet abundance in current scope
        new_triplets: triplet abundance in new scope
    Returns:
        spectrum inferred from mut_catalogue as seen in new scope
    """

    new_profile = {sample: None for sample in mut_catalogue}
    for sample in mut_catalogue:
        profile = mut_catalogue[sample]
        total_mutations = sum(profile.values())
        norm_profile = utils.normalize_profile(profile, curr_triplets)
        new_scope_profile = utils.denorm_profile(norm_profile, new_triplets)
        new_profile[sample] = {c: abs(new_scope_profile[c] * total_mutations) for c in new_scope_profile}
    return new_profile


def prepare_sigfit(annotmuts, gw_triplets, exon_triplets, mut_catalogue_path):
    """
    Args:
        annotmuts: dndsout$annotmuts produced by dNdScv
        gw_triplets: triplet abundance in genome scope
        exon_triplets: triplet abundance in exome scope
        mut_catalogue_path: path where storing mutational catalogue
    Returns:
         mutational catalogue in sigfit format
         cosmic signatures into exome scope
    """
    reader = Reader()
    mut_catalogue = create_mut_catalogue(annotmuts)
    mut_catalogue.to_csv(mut_catalogue_path, sep='\t')
    cosmic_signatures = reader.read_cosmic()
    cosmic_exome = cosmic_sigfit_format(switch_scope(cosmic_signatures, gw_triplets, exon_triplets))
    cosmic_exome.to_csv(reader.cosmic_exome_path, sep='\t')


def cosmic_norm_combine(arr, cosmic, gw_triplets):
    """
    Args:
        arr: array of coefficients (exposures)
        cosmic: cosmic signatures
        gw_triplets: genomewide triplet abundance
    Returns:
        normalized linear combination of signatures
    """
    assert(len(arr) == 30)
    combined_profile = {mk: 0 for mk in utils.mut_key_generator()}
    for i, signature_name in enumerate(cosmic):
        profile = cosmic[signature_name]
        for mk in utils.mut_key_generator():
            combined_profile[mk] += arr[i] * profile[mk]
    total = sum(combined_profile.values())
    combined_profile = {k: combined_profile[k]/total for k in combined_profile}
    norm_comb_profile = utils.normalize_profile(combined_profile, gw_triplets)
    return norm_comb_profile


def spread_expected(sample_expect, contrib_dict):
    """
    Args:
        sample_expect: dict: sample -> expected number of syn mutations
        contrib_dict: dict: sample -> profile in local synonymous scope
    Return:
        expected_dict: dict: sample -> dict(c -> expected number of mutations in context c)
    """
    expected_dict = {}
    #  sample_expect has been derived from the raw mutations file,
    #  so it may include some samples which only harbour mutations in non-SNV elements,
    #  the probabilities of which do not matter to us at the moment
    for sample, expect in sample_expect.items():
        if sample in contrib_dict:
            expected_dict[sample] = {c: contrib_dict[sample][c] * expect for c in contrib_dict[sample]}
    return expected_dict


def gene_predict(genemuts, ensembl_gene, ensembl_hugo_dict):
    """
    Args:
        genemuts: dndsout$genemuts
        ensembl_gene: str: ENSEMBL ID
        ensembl_hugo_dict: dict: ENSEMBL -> HUGO
    :return:
        predicted syn mutations predicted by NB regression in gene with ENSEMBL ID at cohort level
    """
    return genemuts[genemuts['gene_name'] == ensembl_hugo_dict[ensembl_gene]]['exp_syn_cv'].values[0]


def compute_expected_dict(contrib_dict, ensembl_gene, genemuts, rel_syn_burden, ensembl_hugo_dict):
    gene_predicted = gene_predict(genemuts, ensembl_gene, ensembl_hugo_dict)
    sample_expect = {sample: rel_syn_burden[sample] * gene_predicted for sample in rel_syn_burden}
    return spread_expected(sample_expect, contrib_dict)


def compute_bp_expected(expected_dict, syn_sites):
    """
    Args:
        expected_dict: dict: sample -> dict(c -> expected number of mutations in context c)
        syn_sites: dict: dict(c -> number of syn_sites in context c)
    Return:
        bp_expected: dict(sample -> dict(c -> expected number of mutations per bp in context c)
    """
    bp_expected = {}
    for sample in expected_dict:
        bp_expected[sample] = [expected_dict[sample][c] / syn_sites[c] for c in utils.mut_key_generator()]
    return bp_expected


""" main functions """


def genewise_run(hugo_gene, hugo_ensembl_dict, site_count, cosmic, gw_triplets,
                 exposures, genemuts, rel_syn_burden, ensembl_hugo_dict):
    try:
        ensembl_gene = hugo_ensembl_dict[hugo_gene]
        syn_sites = retrieve_counts(site_count, ensembl_gene)
        total = sum(syn_sites.values())
        subs_abundance = {k: syn_sites[k] / total for k in syn_sites}
        contrib_dict = {ind: utils.denorm_subs(cosmic_norm_combine(row.values, cosmic, gw_triplets), subs_abundance)
                        for ind, row in exposures.iterrows()}
        expected_dict = compute_expected_dict(contrib_dict, ensembl_gene, genemuts, rel_syn_burden, ensembl_hugo_dict)
        return {hugo_gene: compute_bp_expected(expected_dict, syn_sites)}
    except:
        #  KeyErrors are expected where an ENSEMBL symbol does not match any keys in the ENSEMBL -> HUGO dictionary
        #  bare except aims to prevent similar errors that may or may not have showed up during tests
        #  like ValueErrors or IndexErrors
        return {hugo_gene: None}


# cli options


@click.group()
def mutrate():
    pass


@mutrate.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--annotmuts_path', '-a',
              help='path for dndscv output: dndsout$annotmuts',
              type=click.Path())
@click.option('--output', '-o', 'mut_catalogue_path',
              help='path for output mutational catalogue',
              type=click.Path())
def sigfit(annotmuts_path, mut_catalogue_path):
    reader = Reader()
    annotmuts = pd.read_csv(annotmuts_path, sep='\t', dtype={'sampleID': 'object'})
    exon_triplets, genomewide_triplets = tuple(map(reader.read_triplets, ['exon', 'full']))
    prepare_sigfit(annotmuts, genomewide_triplets, exon_triplets, mut_catalogue_path)


@mutrate.command(context_settings={'help_option_names': ['-h', '--help']})
@click.option('--annotmuts', '-a', 'annotmuts_path',
              help='path for dndscv output: dndsout$annotmuts',
              type=click.Path())
@click.option('--genemuts', '-g', 'genemuts_path',
              help='path for dndscv output: dndsout$genemuts',
              type=click.Path())
@click.option('--mutcat', '-m', 'mut_catalogue_path',
              help='path for mutational catalogue',
              type=click.Path())
@click.option('--exposures', '-e', 'sigfit_exposures_path',
              help='path for COSMIC exposures coming after sigfit',
              type=click.Path())
@click.option('--cores', '-c', default=os.cpu_count(),
              help='Maximum processes to run in parallel',
              type=click.INT)
@click.option('--output_path', '-o', 'output_path', required=True,
              help='Path to folder where keeping json dictionaries: {gene : {sample : {context : expectation}}}',
              type=click.Path())
def compute_mutrate(annotmuts_path, genemuts_path, mut_catalogue_path, sigfit_exposures_path, cores, output_path):

    # instantiate reader
    reader = Reader()

    # retrieve annotations
    site_count = reader.read_site_counts()
    cosmic = reader.read_cosmic()
    annotmuts = pd.read_csv(annotmuts_path, sep='\t')
    genemuts = pd.read_csv(genemuts_path, sep='\t')
    exon_triplets, gw_triplets = tuple(map(reader.read_triplets, ['exon', 'full']))
    exposures = reader.read_exposures(mut_catalogue_path, sigfit_exposures_path)
    hugo_ensembl_dict, ensembl_hugo_dict = reader.read_gene_dicts()
    rel_syn_burden = compute_rel_syn_burden(annotmuts)

    # Define the geneset: minimal set required to annotate MAF with expected mutations per bp
    dg = annotmuts[annotmuts['mut'].isin(list('ACGT'))]
    gene_set = dg['gene'].unique()
    # gene_set = ['TP53']

    # instantiate a partial of genewise_run, which encapsulates the main task
    run_task = partial(genewise_run, hugo_ensembl_dict=hugo_ensembl_dict,
                       site_count=site_count, cosmic=cosmic, gw_triplets=gw_triplets,
                       exposures=exposures, genemuts=genemuts, rel_syn_burden=rel_syn_burden,
                       ensembl_hugo_dict=ensembl_hugo_dict)

    # loop task through geneset
    result = {}
    with Pool(cores) as pool:
        for a in tqdm(pool.imap(run_task, gene_set), total=len(gene_set)):
            # option 1: dump to separate json files, one per each gene
            with open(os.path.join(output_path, '{0}.out.json'.format(next(iter(a.keys())))), 'wt') as f_output:
                json.dump(a, f_output)

    # option 2: pickle the results in a single dictionary
    # with gzip.open(output_path, 'wb') as f_output:
    #     pickle.dump(result, f_output)


cli = click.CommandCollection(sources=[mutrate])

if __name__ == '__main__':

    """
    Pipeline
    %%%%%%%%

    Requires:
    ========
    $ source activate signatures_env

    step 1:
    ======
    cli example: python mutrate.py sigfit -a <annotmuts_path>
                                          -o <mut_catalogue_path>
    step 2:
    ======
    cli example: Rscript --vanilla cosmic_exome_fit.r -m <mut_catalogue_path>
                                                      -c <cosmic_signatures_path>
                                                      -o output
    step 3:
    ======
    cli example: python mutrate.py compute_mutrate -a <annotmuts_path>
                                                   -g <genemuts_path>
                                                   -m <mut_catalogue_path>
                                                   -e <sigfit_exposures_path>
                                                   -c 30
                                                   -o <folder>
                                                   
    """
    cli()

