import pandas as pd
import os
import argparse
import subprocess

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Given a MAF of mutations it computes deconstructSigs R package. Then with the output results of '
                    'deconstructSigs R package calculates the probability for each mutation (type) generated by each '
                    'signature using W and H matrices.')

    parser.add_argument('-i', '--input_file',
                        dest="input",
                        action="store",
                        default=None,
                        required=True,
                        help="Input must be a string written in command line specifying the absolute path to the "
                             "mutation list of the cohort")
    parser.add_argument('-o', '--output_path',
                        dest="output",
                        action="store",
                        default=None,
                        required=True,
                        help="Output must be a string written in command line specifying the absolute path for "
                             "the output of the cohort")
    parser.add_argument('-s', '--sign',
                        dest="sign",
                        action="store",
                        default=None,
                        required=False,
                        help="Path to the signature's matrices (cosmic signatures)")

    args = parser.parse_args()

    in_path = args.input
    out_path = os.path.abspath(args.output)
    path_matrix = os.path.join(os.environ['INTOGEN_DATASETS'], 'deconstructsig', 'signatures_exons_normalized.cosmic') if args.sign is None else args.sign
    cohort = in_path.split("/")[-1].split(".")[0]

    # READ DATASET AND FILTER OUT INDELS

    muts = pd.read_table(in_path, sep='\t', compression="gzip")
    muts = muts[(muts['REF'] != '-') & (muts['REF'] == 'A') | (muts['REF'] == 'C') | (muts['REF'] == 'T') | (
            muts['REF'] == 'G')]
    muts = muts[(muts['ALT'] != '-') & (muts['ALT'] == 'A') | (muts['ALT'] == 'C') | (muts['ALT'] == 'T') | (
            muts['ALT'] == 'G')]

    # READ DRIVERS FOR THIS COHORT
    DRIVERS = os.path.join(os.environ['INTOGEN_DATASETS'], 'deconstructsig', 'output_pass_drivers_01.csv')
    df_filtered = pd.read_csv(DRIVERS, sep="\t")
    df_filtered = df_filtered[df_filtered["COHORT"] == cohort]

    # FILTER OUT MUTATIONS IN DRIVERS, IF THERE IS ANY
    if df_filtered.shape[0] > 0:
        list_drivers = df_filtered["SYMBOL"].unique()
        muts = muts[~(muts["GENE"].isin(list_drivers))]

    # CREATE OUTPUT FOLDER
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # CREATE TEMPORAL FILE OF MUTATIONS WITHOUT INDELS
    if not (os.path.exists(os.path.join(out_path, "tmp"))):
        os.mkdir(os.path.join(out_path, "tmp"))

    muts.to_csv(os.path.join(out_path, "tmp", "mutations_snv.csv"), sep='\t', index=False)

    # RUN DECONSTRUCTSIGN IN R

    command = 'cd ' + out_path + ' && Rscript ' + os.path.join(os.path.abspath(os.path.dirname(__file__)), "deconstructSigs.r")+ ' ' + os.path.join(out_path, "tmp", "mutations_snv.csv")
    print(command)
    subprocess.run(command, shell=True, executable='/bin/bash')

    # REMOVE TEMPORAL FILE

    os.remove(os.path.join(out_path, "tmp", "mutations_snv.csv"))
    os.rmdir(os.path.join(out_path, "tmp"))

    # READ SIGNATURE FILE (WEIGHT MATRIX)
    W = pd.read_csv(path_matrix, header=0, sep="\t", index_col=0)

    #  TRANSPOSE
    W = W.T

    # READ EXPOSURE FILE
    H = pd.read_csv(os.path.join(out_path, "signatures_weight.csv"), header=0, sep="\t")

    # COMPUTE PROBABILITIES

    # go over each sample in H matrix and compute the probability for each mutation type in a tri-nuc context

    frames = []  # to collect results sample wise
    flag = 0
    for idx, row in H.iterrows():  # go over each sample
        sample = row['sample_id']
        sig_dic = {}
        allsigs = []

        # get the exposure (i.e total number of mutations belong to each signature) value for the particular sample from H matrix
        for col in H.columns:
            if col not in ['sample_id', 'SSE', 'mutation_count']:
                sig_dic[col] = row[col] * row[
                    'mutation_count']  # save the exposuse value in a dictionary per signature name
                allsigs.append(col)  # save the signature names

        # multiple the exposure (from H) with the W matrix
        a = W.copy()  # take a copy of the W matrix (which is the extracted signatures - not sample specific)
        for sig in allsigs:
            a[sig] *= sig_dic[
                sig]  # mutiply the signature columns with the corresponding signature exposure in that particular sample

        # compute the row sum for normalization (i.e sum of values across signature for each mutation/context type)
        a['row_sum'] = a[allsigs].sum(axis=1)

        # normalize the row values with the row sum to driver
        # the probabilities for different signatures for each mutation type
        new = a[allsigs].div(a['row_sum'], axis=0)[allsigs]

        # add info columns
        new['Mutation_type'] = new.index
        new['Sample'] = sample

        # sort the columns
        columns = ['Sample', 'Mutation_type'] + allsigs

        new = new[columns]

        # save the results for each samples in a dataframe
        if flag == 0:
            frames = [new]
            flag += 1
        else:
            frames.append(new)

    results_new = pd.concat(frames)

    # WRITE RESULTS

    # this output file will contains the probabilities for each mutation type (under particular tri-nucleotides) to be
    # generated by different signatures per sample.
    results_new.to_csv(os.path.join(out_path, "mutation_sign_prob.tsv"), sep="\t", header=True, index=False)
