import os

import pandas as pd


def read_file(filein):
    f = open(filein, 'r')
    genes = set()
    for line in f.readlines():
        line = line.strip()
        genes.add(line)
    f.close()
    return genes


def check_black_white_lists(df_all):

    # Read vetted file
    df_drivers = df_all[df_all["FILTER"] == "PASS"]
    print("Number of drivers before white/black listing:" + str(len(df_drivers["SYMBOL"].unique())))

    # Remove black listed genes
    black_listed_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'postprocess', 'black_listed.txt')
    black_listed = read_file(black_listed_file)
    df_drivers = df_drivers[~df_drivers["SYMBOL"].isin(black_listed)]

    # Now rescue white listed discarded genes
    df_discarded = df_all[df_all["FILTER"] == "Lack of literature evidence"]

    # Recover white listed genes
    white_listed_file = os.path.join(os.environ['INTOGEN_DATASETS'], 'postprocess', 'white_listed.txt')
    white_listed = read_file(white_listed_file)
    df_recovered = df_discarded[df_discarded["SYMBOL"].isin(white_listed)]
    df_recovered["FILTER"] = "PASS"

    df_final_list = pd.concat([df_drivers, df_recovered], sort=True)
    print("Number of drivers after white/black listing:" + str(
        len(df_final_list[df_final_list["FILTER"] == "PASS"]["SYMBOL"].unique())))

    return df_final_list