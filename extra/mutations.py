
import pandas as pd


def get_protein_mutation(pos, change):
    wt_aa = change[0]
    mt_aa = change[-1]
    if mt_aa != wt_aa:
        return wt_aa + str(pos) + mt_aa
    else:
        return wt_aa + str(pos)


def information(mutations, driver_genes):
    df = pd.read_csv(mutations, sep="\t")
    df["MUTATION"] = df['Location'].astype(str) + ':' + df['ref'] + '>' + df['alt']
    df["CONSEQUENCE_SELECTED"] = df['Consequence'].str.split(',')[0]
    df['PROTEIN_MUTATION'] = df.apply(
        lambda row: get_protein_mutation(row["Protein_position"], row["Amino_acids"]), axis=1)
    df_g = df.groupby(
        ["MUTATION", "ENSEMBL_TRANSCRIPT", "gene", "COHORT", "CONSEQUENCE_SELECTED", "chr", "pos", "alt",
         "Protein_position", "PROTEIN_MUTATION"], as_index=False).agg({"sample_id": "count"})
    df_g.rename(columns={"sample_id": "SAMPLES"}, inplace=True)
    df_g[["gene", "MUTATION", "ENSEMBL_TRANSCRIPT", "CONSEQUENCE_SELECTED", "Protein_position",
              "Protein_mutation"]].drop_duplicates().to_csv(
        "/workspace/projects/intogen_2017/test/boostDM/web_data/info_mutations_all.tsv.gz", compression="gzip",
        sep="\t", index=False)