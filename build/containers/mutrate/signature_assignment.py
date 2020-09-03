
import click
import pandas as pd

from config import Values

# global variables
COSMIC_EXOME = Values.COSMIC_EXOME


# run functions

@click.command()
@click.option('--input_file', type=click.Path())  # weights from deconstructSigs
# @click.option('--cosmic_exome', type=click.Path())
@click.option('--assignment_path', type=click.Path())
def run_assign(input_file, assignment_path):

    """
    Takes the results of deconstructSigs.r and computes the likelihood that
    a mutation for each context occurs given the sample and the signature.
    """

    # read COSMIC signature file and transpose
    W = pd.read_csv(COSMIC_EXOME, sep="\t", index_col=0, header=0)
    W = W.T

    # read weights
    H = pd.read_csv(input_file, header=0, sep='\t')

    # compute probabilities
    # go over each sample in H matrix and compute the probability for each mutation type in a tri-nuc context

    frames = []  # to collect results sample-wise
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
    results_new.to_csv(assignment_path, sep='\t')


if __name__ == '__main__':
    run_assign()
