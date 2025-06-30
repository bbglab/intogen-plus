'''
Created on May 3, 2017

@author: root
'''
"""
Script to enable multiple testing correction
for hotspots

authors : Collin Tokheim, Rohit Bhattacharya
emails : collintokheim@gmail.com, rohit.bhattachar@gmail.com
"""
# imports
import src.statistics as mystats
import argparse
import operator
import get_hotspot_residues as get_hotspot
import numpy as np
import csv
import os

MISSING_COUNT = 0
EMPTY_PVALS = 0

def parse_arguments():
    """
    Function to parse command line arguements
    from the user

    Returns
    -------
    opts : dict
        command line arguements from the user
    """

    info = 'Multiple testing correction for predicted hotspots'
    parser = argparse.ArgumentParser(description=info)

    # program arguments
    parser.add_argument('-i', '--hotspot-file',
                        type=str,
                        required=True,
                        help='File containing the verbose output from hotspot.py')
    parser.add_argument('-f', '--function',
                        type=str,
                        default='min',
                        help='Function applied to a group of p-values (Default: min)')
    parser.add_argument('-m', '--mupit-dir',
                        type=str,
                        required=True,
                        help='Directory containing the mupit annotations')
    parser.add_argument('-q', '--q-value',
                        type=float,
                        required=True,
                        help='Q-value (FDR) threshold to apply for significance')
    parser.add_argument('-o', '--output-file',
                        type=str,
                        help='Name of output file')
    parser.add_argument('-s', '--significance-level',
                        type=str, required=True,
                        help='Output file to write the significance level cutoff'
                        ' for each tumor type')

    args = parser.parse_args()
    opts = vars(args)
    return opts


def read_mupit_file(in_file):
    """
    Reads in a mupit annotated file and sorts it
    by gene, transcript and residue

    Parameters:
    -----------
        in_file : str
            name of the mupit annotated file

    Returns:
    --------
        mupit_annotations : list
            sorted by gene, transcript and residue
    """

    # mupit annotations list
    mupit_annotations = []

    # open the file
    in_file = open(in_file)

    # read in header
    header = in_file.readline().strip()
    header = header.split('\t')

    # get indices of relevant info
    gene_ind = header.index("HUGO symbol")
    #transcript_ind = header.index("Sequence ontology transcript")
    #residue_ind = header.index("Sequence ontology protein sequence change")
    transcript_ind = header.index("Reference Transcript")
    residue_ind = header.index("Reference Codon Position")
    aa_ind = header.index("Reference AA")
    genomic_pos_ind = header.index("Reference Genomic Position")
    chrom_ind = header.index("Chromosome")

    for line in in_file:
        # remove trailing characters
        line = line.strip()
        # split by tab
        line = line.split('\t')

        # obtain the cravat residue id, which is between two letters indicating
        # the protein sequence change
        # also leave out special cases
        if len(line) < residue_ind or line[residue_ind] == '':
            continue

        # get residue position
        #line[residue_ind] = int(line[residue_ind][1:-1])
        line[residue_ind] = int(line[residue_ind])

        # add as tuple so we can sort
        mupit_annotations.append(tuple(line))

    # sort by gene, transcript and position
    mupit_annotations = sorted(mupit_annotations,
                               key=operator.itemgetter(gene_ind, transcript_ind, residue_ind))
    in_file.close()
    return (mupit_annotations, gene_ind, transcript_ind,
            residue_ind, chrom_ind, genomic_pos_ind, aa_ind)


def get_group_pvals(mupit_groups, gene_ind, transcript_ind,
                    residue_ind, chrom_ind, genomic_pos_ind, aa_ind,
                    hotspot_output, ttype, func=min):
    """
    """

    global MISSING_COUNT
    global EMPTY_PVALS


    # dictionary containing pdb id, chain, residue as keys
    # and a list of all pvals associated with the unique combination
    # as values
    hspot_pvals = {}

    # iterate through each line of the hotspot output for the given ttype
    for hspot_data in hotspot_output:
        # make a key from pdb id, chain, and residue
        curr_key = (hspot_data[0], hspot_data[3], hspot_data[4])

        # check if this key is already in our pvals dict
        # if not, add it
        if not curr_key in hspot_pvals:
            hspot_pvals[curr_key] = [float(hspot_data[5])]
        else:
            hspot_pvals[curr_key].append(float(hspot_data[5]))

    # list of pvals grouped by unique gene, transcript, residue
    grouped_pvals = []

    # list of min pvals from each grouping of pvals
    min_pvals = []

    # variable to check last grouping of gene, transcript, residue
    prev_group = (mupit_groups[0][gene_ind],
                  mupit_groups[0][transcript_ind],
                  mupit_groups[0][residue_ind])
    prev_line = mupit_groups[0]

    # go through all mupit data
    for j, line in enumerate(mupit_groups):

        # get current group
        curr_group = (line[gene_ind], line[transcript_ind], line[residue_ind])

        # check if this is the same as previous
        if not prev_group[0] == curr_group[0] \
           or not prev_group[1] == curr_group[1] \
           or not prev_group[2] == curr_group[2]:

            # check if we have groups with zero
            # pvals
            if not grouped_pvals:
                EMPTY_PVALS += 1
                #print "no pval"
            else:
                # if not, get the min pval and move onto the next
                min_pvals.append([prev_group[0], ttype, prev_group[1],
                                  prev_group[2], prev_line[aa_ind], prev_line[chrom_ind],
                                  prev_line[genomic_pos_ind], func(grouped_pvals)])
                grouped_pvals = []

        # update prev_group
        prev_group = curr_group
        # update prev_line
        prev_line = line

        # get the pvals corresponding to the pdb, chain, residue
        # combination, stored in hspot_pvals dict
        pdb_chain_residue = (line[0], line[1], line[2])

        # check if this combination exists in our hotspot data
        if not pdb_chain_residue in hspot_pvals:
            MISSING_COUNT += 1
            continue

        # add all associated pvals
        for pval in hspot_pvals[pdb_chain_residue]:
            grouped_pvals.append(pval)

    return min_pvals


def main(opts):
    """
    Main function
    """
    # obtain user defined parameters
    hspots_file = opts["hotspot_file"]
    out_file = opts["output_file"]
    mupit_dir = opts["mupit_dir"]
    signif_lvl_file = opts['significance_level']

    # use external module to separate out the residues in the hotspot.py output
    # onto separate lines
    args = {"input": hspots_file,
            "significance_level": 1.1,
            "output": None}
    hotspot_output = get_hotspot.main(args)

    # stratify hotspot output by tumour type
    stratified_hotspot_output = {}
    for hspot_data in hotspot_output:
        # skip homology model
        #if hspot_data[0].startswith('NP_') or hspot_data[0].startswith('ENSP'):
            #continue

        ttype = hspot_data[1]
        if not ttype in stratified_hotspot_output:
            stratified_hotspot_output[ttype] = [hspot_data]
        else:
            stratified_hotspot_output[ttype].append(hspot_data)

    # open a file to write results to
    header = ["HUGO Symbol", "Tumor Type",
              "Sequence Ontology Transcript", "CRAVAT Res",
              "Ref AA",
              'chromosome', 'genomic position',
              "Min p-value", 'q-value']
    output, signif_lvl_output = [], []
    # go through all mupit annotation files
    if opts['function'] == 'min':
        myfunc = min
    elif opts['function'] == 'median':
        myfunc = np.median
    elif opts['function'] == 'max':
        myfunc = max

    # HACK for m_file in os.listdir(mupit_dir):
    m_file = mupit_dir
    # print the tumor type we're on
    # ttype = m_file.split('_')[-1]
    ttype = m_file.split('/')[-1].split(".")[0]
    # skip cancer type if no mutations/hotspots.
    # usually only happens if they use only a couple structures
    # and not the full PDB
    if ttype in stratified_hotspot_output:
        print(ttype)

        # read in the file
        tmp = read_mupit_file(m_file)
        (mupit_annotations, gene_ind, transcript_ind,
         residue_ind, chrom_ind, genomic_pos_ind, aa_ind) = tmp

        # get p-values for grouped up mupit-hotspot groups
        grouped_p_vals = get_group_pvals(mupit_annotations, gene_ind,
                                         transcript_ind, residue_ind,
                                         chrom_ind, genomic_pos_ind, aa_ind,
                                         stratified_hotspot_output[ttype],
                                         ttype, myfunc)
        # add the q-value
        tmp_pvals = [g[-1] for g in grouped_p_vals]
        tmp_qvals = mystats.bh_fdr(tmp_pvals)
        for i in range(len(grouped_p_vals)):
            grouped_p_vals[i].append(tmp_qvals[i])
        output.extend(grouped_p_vals)

        # figure out the equivalent p-value threshold
        signif_pvals = [tmp_pvals[i]
                        for i in range(len(tmp_qvals))
                        if tmp_qvals[i] <= opts['q_value']]
        if signif_pvals:
            pval_cutoff = max(signif_pvals)
            signif_lvl_output.append([ttype, pval_cutoff])
    # HACK: END FOR

    # write to output file
    with open(out_file, 'wb') as out_file:
        mywriter = csv.writer(out_file, delimiter='\t', lineterminator='\n')
        mywriter.writerows([header]+output)

    # write significance level cutoffs for each tumor type
    with open(signif_lvl_file, 'wb') as out_file:
        mywriter = csv.writer(out_file, delimiter='\t', lineterminator='\n')
        mywriter.writerows(signif_lvl_output)
    #print("MISSING COUNTS = " + str(MISSING_COUNT))


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
