import click
import pandas as pd
import numpy as np

def concat(grp):
    return ",".join([str(x) for x in list(grp)])

def checkMNVS(li):
    """ Check if any of the position in the list are consecutive"""
    l = [int(v) for v in li]
    n = len(l) - 1
    return (sum(np.diff(sorted(l)) == 1) >= n)

def run(vep_files):
    list_dfs = []
    for filein in vep_files:
        cohort = filein.split(".")[0]
        try:
            df = pd.read_csv(filein, sep="\t", low_memory=False)
        except:
            continue
        df["COHORT"] = cohort
        ##Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACTDISTANCE	STRAND	FLAGS	SYMBOL	SYMBOL_SOURCE	HGNC_ID	CANONICAL	ENSP
        df = df[(df["CANONICAL"]=="YES")][["#Uploaded_variation","Location","Feature","SYMBOL","COHORT"]]
        
        df[["id_", "sample_id", "ref", "alt", "pos"]] = df["#Uploaded_variation"].str.split("__", expand=True)
        df["chr"] = df["Location"].str.split(":", expand=True)[1]

        df=df.groupby(["SYMBOL","sample_id","COHORT"],as_index=False).agg({"pos":concat})
        
        df["MNV"] = df[df["pos"].str.contains(",")].pos.str.split(',').apply(lambda x: checkMNVS(x))
        df["MNV"].fillna(False, inplace=True)
        
        list_dfs.append(df.drop_duplicates().copy())
    
    df_fin = pd.concat(list_dfs)

    df_fin.to_csv("mnvs.tsv.gz", sep="\t", index=False, compression="gzip")


@click.command(context_settings={'help_option_names': ['-h', '--help']})
@click.argument('vep_files', nargs=-1)
def cli(vep_files):
    run(vep_files)
    

if __name__ == "__main__":
    cli()