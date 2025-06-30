import csv
import numpy as np
import re
import time
from collections import Counter

def calculate_cutoff(samp_list, mut_threshold=None):
    """Calculates Kandoth et al. hypermutator cutoff or just a single
    mutation count definition for hypermutator status."""
    mut_cts = Counter(samp_list)
    mut_vals = mut_cts.values()
    if mut_threshold is None:
        q75, q25 = np.percentile(mut_vals, [75 ,25])
        iqr = q75 - q25
        mut_cutoff = q75 + iqr*4.5
    else:
        mut_cutoff = mut_threshold
    return mut_cts, mut_cutoff


def detect_hypermutators(file_path,
                         samp_colname='Tumor_Sample_Barcode',
                         stratify_col=None,
                         mut_threshold=None):
    """Finds which samples are hyper-mutators based on 4.5 * IQR.

    Threshold based on Kandoth et al. paper.

    Parameters
    ----------
    file_path : str
        path to TCGA mutation file
    samp_colname : str
        column name with sample IDs
    stratify_col : str or None
        column for stratified hyper-mutator detection. Likely indicating
        tumor types.
    mut_threshold : int or None
        A flat number of mutations to filter hypermutators. This will
        overide the quantile based evaluation.

    Returns
    -------
    hypermut_samps : list
        sample IDs that are hypermutators
    num_mut_list : list
        number of mutations for each hypermutator sample
    """
    # get mutation counts for each sample
    with open(file_path) as handle:
        myreader = csv.reader(handle, delimiter='\t')

        # parse header
        while True:
            header = next(myreader)
            if header[0].startswith('#'):
                # skip comment lines
                continue
            else:
                samp_ix = header.index(samp_colname)
                break

        # get hypermutator cutoff based on Kandoth et al.
        if not stratify_col:
            data = list(myreader)
            samp_list = [t[samp_ix] for t in data]
            mut_cts, mut_cutoff = calculate_cutoff(samp_list, mut_threshold)
        else:
            data = list(myreader)
            samp2grp = dict()
            strat_ix = header.index(stratify_col)

            # figure out what strat each sample corresponds to
            for line in data:
                samp2grp[line[samp_ix]] = line[strat_ix]

            # get hyper-mutator threshold for each strata
            grps = set(t[strat_ix] for t in data)
            grp2cutoff = dict()
            for grp in grps:
                tmp_data = [t for t in data if t[strat_ix]==grp]
                samp_list = [t[samp_ix] for t in tmp_data]
                if not samp_list:
                    #IPython.embed()
                    pass
                mut_cts, mut_cutoff = calculate_cutoff(samp_list, mut_threshold)
                grp2cutoff[grp] = mut_cutoff

        # classify hypermutator samples
        hypermut_samps = ['hypermutated sample id']
        num_mut_list = ['# mutations']
        mut_cts = Counter([d[samp_ix] for d in data])
        for sample in mut_cts:
            tmp_num_mut = mut_cts[sample]
            if not stratify_col:
                if tmp_num_mut > mut_cutoff:
                    hypermut_samps.append(sample)
                    num_mut_list.append(tmp_num_mut)
            else:
                tmp_mut_cutoff = grp2cutoff[samp2grp[sample]]
                if tmp_num_mut > mut_cutoff:
                    hypermut_samps.append(sample)
                    num_mut_list.append(tmp_num_mut)

    return hypermut_samps, num_mut_list


def comment_stripper(iterator):
    """Skips over commented or blank lines."""
    for line in iterator:
        if line [:1] == '#':
            continue
        if not line.strip():
            continue
        yield line
