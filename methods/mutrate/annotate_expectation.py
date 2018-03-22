import os
import pandas as pd
import itertools
import json
import click
import glob
import gzip
from contextlib import contextmanager
from tqdm import tqdm


CB = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


@contextmanager
def generic_open_file(file):
    if file.endswith('.gz'):
        f = gzip.open(file, 'rt')
    else:
        f = open(file, 'rt')
    yield f
    f.close()

def mut_key_generator():

    subs = ['CA', 'CG', 'CT', 'TA', 'TC', 'TG']
    for s in sorted(subs):
        for c in sorted(itertools.product({'A', 'C', 'G', 'T'}, repeat=2)):
            yield tuple([s, ''.join(c)])


def purine_to_pyrimidine(m):
    if m[0][0] not in list('CT'):
        t1 = CB[m[0][0]] + CB[m[0][1]]
        t2 = CB[m[1][1]] + CB[m[1][0]]
        return tuple([t1, t2])
    else:
        return m


mut_keys = [k for k in mut_key_generator()]


def retrieve_expectation(exp_dict, v):
    sample = v['sampleID']
    ref3_cod = v['ref3_cod']
    mut3_cod = v['mut3_cod']
    key = tuple([ref3_cod[1] + mut3_cod[1], ref3_cod[0] + ref3_cod[2]])
    key = purine_to_pyrimidine(key)
    return exp_dict[sample][mut_keys.index(key)]


@click.command(context_settings={'help_option_names': ['-h', '--help'], 'max_content_width': 120})
@click.option('--genes_path', type=click.Path(exists=True), help='folder with gene jsons')
@click.option('--maf_path', type=click.Path(exists=True), help='mutations file')
@click.option('--dnds_path', type=click.Path(exists=True), help='dnds output file')
@click.option('--output', '-o', type=click.Path(exists=False), help='output_file')
def mutrate(genes_path, maf_path, dnds_path, output):
    """
    \b
    python annotate_expectation.py \\
    --genes_path /workspace/projects/intogen_2017/runs/20180307/mutrate/PCATLAS_WXS_UCEC/genes/ \\
    --maf_path /workspace/projects/intogen_2017/runs/20180307/dndscv/PCATLAS_WXS_UCEC_annotmuts.out.gz \\
    --dnds_path /workspace/projects/intogen_2017/runs/20180307/dndscv/PCATLAS_WXS_UCEC.out.gz.tsv \\
    --output /workspace/users/fmuinos/mutrate/ppln_results/BRCA_TCGA/UCEC_annotmuts_expect.tsv
    """
    annotmuts = pd.read_csv(maf_path, sep='\t')
    dnds = pd.read_csv(dnds_path, sep='\t')
    assert('mut' in annotmuts.columns)
    assert('ref' in annotmuts.columns)
    assert('gene' in annotmuts.columns)
    annotmuts = annotmuts[(annotmuts['ref'].isin(list('ACGT'))) & (annotmuts['mut'].isin(list('ACGT')))]
    df = pd.DataFrame({})
    for gene in tqdm(annotmuts['gene'].unique()):
        gene_file_path = glob.glob('{0}/{1}.*'.format(genes_path, gene)).pop()
        try:
            with generic_open_file(gene_file_path) as f:
                gene_json_dict = json.load(f)
            dg = annotmuts[annotmuts['gene'] == gene]
            dg['expect'] = dg.apply(lambda v: retrieve_expectation(gene_json_dict[gene], v), axis=1)
            df = pd.concat([df, dg])
        except:
            pass
    merged = df.merge(dnds[['gene_name', 'wmis_cv', 'wnon_cv', 'wspl_cv', 'pallsubs_cv', 'qallsubs_cv']], left_on='gene', right_on='gene_name', how='outer').dropna()
    merged.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    mutrate()


