"""Filters out hyper-mutated samples."""
import os
from collections import Counter
import numpy as np
import csv
import re
import argparse


def parse_arguments():
    info = 'Filters out hypermutated samples'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-r', '--raw-dir',
                        type=str, default='/home/pipeline/mupit_update/tcga/',
                        help='Directory where Xena data/Raw data is located')
    parser.add_argument('-m', '--match-regex',
                        type=str, default='^TCGA.+Xena\.dat$',
                        help='Regular expression matching file name of interests in the raw directory')
    parser.add_argument('-s', '--sample-col',
                        type=str, default='#sample',
                        help='Column name for sample ID')
    parser.add_argument('-t', '--tumor-type-col',
                        type=str, default=None,
                        help='Column name indicating tumor type if using PANCAN '
                        'data (optional).')
    parser.add_argument('-mt', '--mut-threshold',
                        type=int, default=None,
                        help='Use flat mutation count for hypermutator filter (Default: None)')
    parser.add_argument('-d', '--data-dir',
                        type=str, default='/home/pipeline/mupit_update/tcga/',
                        help='Directory to save result (Default: /home/pipeline/mupit_update/tcga/)')
    args = parser.parse_args()
    return vars(args)


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
                         samp_colname='#sample',
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


def main(opts):
    # first get hypermutators
    for filename in os.listdir(opts['raw_dir']):
        if opts['match_regex'] == filename:
            mut_file_path = os.path.join(opts['raw_dir'], filename)
            tissue = filename.split('.')[1]
            '''
            if tissue != 'PANCAN12':
                # detect which samples are hyper-mutators
                #hypermut_samps, num_mut_list = detect_hypermutators(mut_file_path,opts['sample_col'],opts['tumor_type_col'],opts['mut_threshold'])
                # write the hypermutated samples to a file
                hypermut_file = os.path.join(opts['data_dir'], 'hypermutated.'+tissue+'.txt')
                with open(hypermut_file, 'w') as writer:
                    for i, l in enumerate(hypermut_samps):
                        writer.write(l+'\t'+str(num_mut_list[i])+'\n')
            else:
                hypermut_samps = []
            '''
            hypermut_samps = []
            # get file path to read in
            #mupit_filename = 'non_filtered_mupit.TCGA.{tissue}.Xena.dat'.format(tissue=tissue)
            mupit_filename = 'non_filtered_mupit.' + filename
            mupit_filepath = os.path.join(opts['data_dir'], mupit_filename)
            if tissue != 'PANCAN12':
                print 'Working on {0}'.format(tissue)
            else:
                print '{0} is already filtered. Will skip filter.'.format(tissue)

            # filter out lines with hypermutated samples
            output = []
            with open(mupit_filepath) as handle:
                myreader = csv.reader(handle, delimiter='\t')

                # iterate through each line
                for line in myreader:
                    additional_info = line[2].split(';')
                    if additional_info[0] not in hypermut_samps:
                        output.append(line)

            # write hypermutated filtered output to mupit.* file
            #out_filename = 'mupit.TCGA.{tissue}.Xena.dat'.format(tissue=tissue)
            out_filename = 'mupit.' + filename
            out_filepath = os.path.join(opts['data_dir'], out_filename)
            with open(out_filepath, 'w') as write_handle:
                mywriter = csv.writer(write_handle, delimiter='\t', lineterminator='\n')
                mywriter.writerows(output)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
