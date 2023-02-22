
import os
import pandas as pd
import click


def get_info(variation,mut):    #O(1)
    #I0000000000__00493087-9d9d-40ca-86d5-936f1b951c93__C__A	1:1787655
    id_,sample_id,ref,alt, pos = variation.split("__")
    chr_,pos_ = mut.split(":")
    output = pd.Series([str(chr_),int(pos),ref,alt,sample_id])
    return output


def concat(grp):
    return ",".join([str(x) for x in list(grp)]) # O(grp_n)
    

def get_mnv(pos):       #O(pos_i)
    if ","  in pos:
        p = [int(v) for v in pos.split(",")]        # O(pos_i)
        return sorted(p) == list(range(min(p), max(p)+1))
    else:
        return False


def main(filein, outfile):
    cohort = os.path.basename(filein).split(".")[0]
    try:
        df = pd.read_csv(filein, sep="\t", low_memory=False)
    except:
        list_dfs = []
        
    df["COHORT"] = cohort
    ##Uploaded_variation	Location	Allele	Gene	Feature	Feature_type	Consequence	cDNA_position	CDS_position	Protein_position	Amino_acids	Codons	Existing_variation	IMPACTDISTANCE	STRAND	FLAGS	SYMBOL	SYMBOL_SOURCE	HGNC_ID	CANONICAL	ENSP
    df = df[(df["CANONICAL"]=="YES")][["#Uploaded_variation","Location","Feature","SYMBOL","COHORT"]]
    df[["chr","pos","ref","alt","sample_id"]] = df.apply(lambda row: get_info(row["#Uploaded_variation"],row["Location"]),axis=1) #O(1) x O(row_j)
    df=df.groupby(["SYMBOL","sample_id","COHORT"],as_index=False).agg({"pos":concat}) #O(grp_n)
    df["MNV"] = df.apply(lambda row: get_mnv(row["pos"]),axis=1) #O(row_j) x O(pos_i)
    
    list_dfs.append(df.drop_duplicates().copy()) 
    


@click.command()
@click.option('-i', '--input', type=click.Path(exists=True), required=True)
@click.option('-o', '--output', type=click.Path(), required=True)
def cli(input, output):
    
    main(input, output)


if __name__ == '__main__':
    cli()