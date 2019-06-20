'''
Compute the driver list from the output of intogen and the vetting information
'''

import pandas as pd
pd.options.mode.chained_assignment = None
import glob
import json
import os
import numpy as np
import json
import click

def read_file(filein):
    f = open(filein,'r')
    genes = set()
    for line in f.readlines():
        line = line.strip()
        genes.add(line)
    f.close()
    return genes



def main(vetting_file,black_listed_file,all_drivers,white_listed_file,ensembl_file,threshold,dir_out):

    # Read vetted file
    df_drivers = pd.read_csv(vetting_file,sep="\t")
    print("Number of drivers before white/black listing:" + str(
        len(df_drivers[df_drivers["FILTER"] == "PASS"]["SYMBOL"].unique())))
    # Read file
    black_listed = read_file(black_listed_file)
    # Remove black listed genes
    df_drivers = df_drivers[~df_drivers["SYMBOL"].isin(black_listed)]
    # Now rescue white listed discarded genes
    df_discarded = pd.read_csv(all_drivers,sep="\t")
    df_discarded = df_discarded[df_discarded["FILTER"]=="Lack of literature evidence"]
    # Recover white listed genes
    white_listed = read_file(white_listed_file)
    df_recovered = df_discarded[df_discarded["SYMBOL"].isin(white_listed)]
    df_recovered["FILTER"] = "PASS"

    df_final_list = pd.concat([df_drivers,df_recovered],sort=True)
    df_final_list.to_csv(
        os.path.join(dir_out, "vetted_drivers" + threshold + ".tsv"), sep="\t", index=False)
    print("Number of drivers after white/black listing:" + str(
        len(df_final_list[df_final_list["FILTER"] == "PASS"]["SYMBOL"].unique())))

    # Create a unique file of drivers
    drivers=df_final_list[df_final_list["FILTER"] == "PASS"]["SYMBOL"].unique()
    # Add the ensembl gene id
    df_ensembl = pd.read_csv(ensembl_file, sep="\t", index_col=False, usecols=[0,1,2,10], names=["ENSEMBL_GENE","SYMBOL","ENSEMBL_PROTEIN","ENSEMBL_TRANSCRIPT"], header=None) # ENSG00000160752	FDPS	ENSP00000349078	1	155312255	155312395	340	480	1260	1	ENST00000356657
    df_drivers_unique= df_ensembl[df_ensembl["SYMBOL"].isin(drivers)].drop_duplicates()
    df_drivers_unique.to_csv(os.path.join(dir_out,"unique_drivers"+threshold+".tsv"), sep="\t",index=False)


@click.command()
@click.option('-i', '--drivers_file', 'drivers_file',type=click.Path(),  help='Path to the raw drivers file',required=True)
@click.option('-a', '--all_drivers_file', 'all_drivers_file',type=click.Path(),  help='Path to the raw genes with all drivers info',required=True)
@click.option('-w', '--white_listed', 'white_listed', type=click.Path(),  help='Path to white listed genes', required=True)
@click.option('-b', '--black_listed', 'black_listed', type=click.Path(),  help="Path to black listed genes", required=True)
@click.option('-o', '--output', 'output', type=click.Path(),  help="Path to the output folder.", required=True)
@click.option('-t', '--threshold', 'threshold',  help="Threshold of the output files, 05 or 01. Default 05", default="05")
@click.option('-n', '--ensembl_file', 'ensembl_file', type=click.Path(),  help="Path to the canonical trancript information file",default="/workspace/projects/intogen_2017/pipeline/datasets/hg38_vep92_develop/shared/cds_biomart.tsv")


def cmdline(drivers_file,black_listed,all_drivers_file,white_listed,ensembl_file,threshold,output):
    main(drivers_file,black_listed,all_drivers_file,white_listed,ensembl_file,threshold,output)



if __name__ == "__main__":
    cmdline()

