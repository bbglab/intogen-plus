import sys
import pandas as pd
import numpy as np


def find_hotspots_gene_cancer(file_hotspots_gene):
    # Structure	Cancer_Type	HotSpot1	HotSpot2...
    f = open(file_hotspots_gene)
    list_hotspots = []

    for line in f:
        line = line.rstrip()
        data = line.split('\t')
        gene = data[0]
        cancer = data[1]
        hotspots = data[2:len(data)]
        number = len(hotspots)
        # Loop over the hotspots, gets the average number of residuse per hotspot and the total number of residues involved in hotspots
        set_h = set()
        sizes_hotspots = []
        for hs in hotspots:
            aas = hs.split(";")
            size = 0
            # ENST00000263967:1047;ENST00000263967:1043	 (Transcript:RES)
            for aa in aas:

                if aa in set_h:
                    size = size + 1
                else:
                    size = size + 1
                    set_h.add(aa)
            sizes_hotspots.append(size)
        list_hotspots.append([gene, cancer, number, len(set_h), np.median(sizes_hotspots)])

    f.close()
    return pd.DataFrame(list_hotspots, columns=["GENE", "Cancer_Type", "TOTAL_HSP", "TOTAL_RES_HSP", "MEDIAN_RES_HSP"])


def get_density_and_pvalue_gene(file_info_density_gene):
    # Structure of the input file
    # Structure       Tumor Type      Model   Chain   Mutation Residues       Residue Mutation Count  Mutation Density        Hotspot P-value
    df_density_aa = pd.read_csv(file_info_density_gene,sep="\t")
    return df_density_aa


def get_pvalue(row,df_density_aa):
    if row["Cancer_Type"] == "REF":
        return 0.0
    return np.min(df_density_aa[(df_density_aa["HUGO Symbol"]==row["GENE"])&(df_density_aa["Tumor Type"]==row["Cancer_Type"])]["Min p-value"].values)


def get_qvalue(row,df_density_aa):
    if row["Cancer_Type"] == "REF":
        return 0.0
    return np.min(df_density_aa[(df_density_aa["HUGO Symbol"]==row["GENE"])&(df_density_aa["Tumor Type"]==row["Cancer_Type"])]["q-value"].values)


def main(file_hotspots_gene, f_info_density, f_output, f_output_clusters):
    df_hotspots = find_hotspots_gene_cancer(file_hotspots_gene)
    df_density_aa = get_density_and_pvalue_gene(f_info_density)
    df_density_aa = df_density_aa.groupby(["HUGO Symbol", "Sequence Ontology Transcript", "Tumor Type"],
                                          as_index=False).agg({"Min p-value": np.min, "q-value": np.min})

    if len(df_hotspots) > 0:
        df_hotspots["Min p-value"] = df_hotspots.apply(lambda row: get_pvalue(row, df_density_aa), axis=1)
        df_hotspots["q-value"] = df_hotspots.apply(lambda row: get_qvalue(row, df_density_aa), axis=1)
        df_hotspots = df_hotspots[df_hotspots["Cancer_Type"] != 'REF']
    else:
        df_hotspots["Min p-value"] = np.nan
        df_hotspots["q-value"] = np.nan

    rows = []
    for index, row in df_density_aa.iterrows():
        if df_hotspots[(df_hotspots["GENE"] == row["HUGO Symbol"]) & (df_hotspots["Cancer_Type"] == row["Tumor Type"])].shape[0] == 0:
            rows.append([row["HUGO Symbol"], row["Tumor Type"], 0, 0, 0, row["Min p-value"], row["q-value"]])

    df_non_significant = pd.DataFrame(rows, columns=df_hotspots.columns.values)

    df_hotspots_all = pd.concat([df_hotspots, df_non_significant])


    df_hotspots_all.to_csv(f_output, sep="\t", index=False, compression="gzip")
    df_density_aa.to_csv(f_output_clusters, sep="\t", index=False, compression="gzip")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
