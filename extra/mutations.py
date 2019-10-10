import glob
import sys
from os import path

import pandas as pd


VALID_CONSEQUENCES = {
        "transcript_ablation", "splice_donor_variant", "splice_acceptor_variant", "stop_gained", "frameshift_variant",
        "stop_lost", "initiator_codon_variant", "transcript_amplification", "inframe_insertion", "inframe_deletion",
        "missense_variant", "splice_region_variant", "incomplete_terminal_codon_variant", "stop_retained_variant",
        "synonymous_variant", "coding_sequence_variant"
    }  # Same than those filtered in intogen pipeline


def valid_csqn(series):
    results = []
    for value in series:
        csqns = value.split(",")
        found = False
        for csqn in csqns:
            if csqn in VALID_CONSEQUENCES:
                results.append(True)
                found = True
                break
        if not found:
            results.append(False)
    return results


def get_info(variation, mut):
    #I0000000000__00493087-9d9d-40ca-86d5-936f1b951c93__C__A	1:1787655
    id_mut, sample_id, ref, alt = variation.split("__")
    chr_, pos = mut.split(":")
    output = pd.Series([str(chr_), int(pos), ref, alt, sample_id])
    return output


def get_protein_mutation(pos, change):
    wt_aa = change[0]
    mt_aa = change[-1]
    if mt_aa != wt_aa:
        return wt_aa + str(pos) + mt_aa
    else:
        return wt_aa + str(pos)


def run(paths, output):
    list_dfs = []
    for path_ in paths:
        for filein in glob.glob(path.join(path_, 'vep', "*.out.gz")):
            cohort = path.basename(filein).split(".")[0]
            try:
                df = pd.read_csv(filein, sep="\t", low_memory=False)
            except:
                continue

            ##Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACTDISTANCE	STRAND	FLAGS	SYMBOL	SYMBOL_SOURCE	HGNC_ID	CANONICAL	ENSP
            df = df[(df["CANONICAL"] == "YES") & (valid_csqn(df["Consequence"]))][["#Uploaded_variation", "Location", "Feature", 'Consequence']]
            df[["CHR", "POS", "REF", "ALT", "SAMPLES"]] = df.apply(lambda row: get_info(row["#Uploaded_variation"], row["Location"]), axis=1)
            df.rename(columns={'Feature': 'TRANSCRIPT'}, inplace=True)
            df["COHORT"] = cohort
            df["MUTATION"] = df['CHR'] + ':' + df['POS'].astype(str) + ':' + df['REF'] + '>' + df['ALT']
            df["CONSEQUENCE"] = df['Consequence'].str.split(',').str[0]
            df.drop(labels=["#Uploaded_variation", "Location", "Consequence"], axis=1)
            df = df.groupby(["MUTATION", "COHORT", "CONSEQUENCE", "CHR", "POS", "REF",
                             "ALT", "TRANSCRIPT"], as_index=False).agg({"SAMPLES": "count"})
            list_dfs.append(df.drop_duplicates())

    df_final = pd.concat(list_dfs)
    df_final.to_csv(output, sep='\t', index=False)


if __name__ == "__main__":
    output = sys.argv[1]
    paths = sys.argv[2:]
    run(paths, output)
