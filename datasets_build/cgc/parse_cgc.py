import os
import pandas as pd
import json
import click
import bglogs

def load_data(path_cgc,path_cohorts,dict_mapping_cgc,dict_mapping_cgc_intogen):
    '''
    Returns the dataframe and dictionaries loaded into python objects
    :param path_cgc:
    :param path_cohorts:
    :param dict_mapping_cgc:
    :param dict_mapping_cgc_intogen:
    :return:
    '''
    # Read data
    cgc_dataset = pd.read_csv(path_cgc, sep=',', header=0)
    intogen_ctypes = pd.read_csv(path_cohorts, sep='\t', header=0)
    mapping = json.load(open(dict_mapping_cgc, "rb"))
    mapping_cgc_intogen = json.load(open(dict_mapping_cgc_intogen, 'rb'))
    return cgc_dataset, intogen_ctypes, mapping, mapping_cgc_intogen

def map_cgc_into_acronym(cgc_dataset,mapping):
    '''
    Read the cgc dataset and include a new column that converts the tumor types into acronyms easily to convert to intogen ttypes
    :param cgc_dataset:
    :return:
    '''
    acronyms = []
    for index, row in cgc_dataset.iterrows():
        # field
        orig_ttypes = str(row["Tumour Types(Somatic)"]).strip()
        cancer_types = []
        if str(orig_ttypes) == "nan":
            acronyms.append("")
            continue
        for orig in orig_ttypes.split(","):
            orig = orig.strip()
            if orig != '':
                try:
                    cancer_types.append(mapping[orig])
                except KeyError:
                    bglogs.warning("Warning, tumor type ",orig,"not found in ",mapping,". Please ensure the tumor type is mapped into an acronym.")
                    continue
        acronyms.append(",".join(cancer_types))
    return acronyms

def map_acronyms_into_intogen(acronyms, mapping_cgc_intogen):
    '''
    Read the cgc dataset and include a new column that converts the acronyms into intogen following the next rules:
    RULE 1 - adenomas mapped to carcinomas of the same cell type (PAAD and PAAC, HA and HC)
    RULE 2 - very specific types/subtypes which will never constitute independent cohort mapped to EXISTING intogen cohorts (BRCAL to BRCA, LIP to S, but BLY not to LY: non-existent)
    RULE 3 - "parent type" to existing cohorts in intogen (NSCLC to LUSC and LUAD, R to RCCC, RCH, RCP)
    :param cgc_dataset
    :return the new column
    '''
    print (acronyms)
    cancer_types = []
    if str(acronyms) == "nan" or acronyms == "":
        return []
    for acronym in acronyms.split(","):
        acronym = acronym.strip()
        if acronym != '' and acronym in mapping_cgc_intogen:
            # The acronym is found, check whether is a leave of a parent node
            ttypes = mapping_cgc_intogen[acronym]
            if len(ttypes) == 1 and ttypes[0] == acronym:
                cancer_types+=[acronym]
            else: # There are several mappings
                for ttype in ttypes:
                    if ttype == acronym: # is a leaf node
                        print (ttype,"ey")
                        cancer_types+= [ttype]
                    else:
                        cancer_types += map_acronyms_into_intogen(ttype,mapping_cgc_intogen)

    return list(set(cancer_types))

@click.command()
@click.option('--path_cgc_original',help= 'path to the downloaded file from cgc', type=click.Path(),required=True)
@click.option('--path_cohorts',help= 'path to the intogen information of the avaliable cohorts', type=click.Path(),required=True)
@click.option('--dict_mapping_cgc',help= 'mapping of text of cgc to acronyms', type=click.Path(),required=True)
@click.option('--dict_mapping_cgc_intogen',help= 'mapping cancer types cgc to intogen', type=click.Path(),required=True)
@click.option('--path_output',help= 'path to the output directory', type=click.Path(),required=True)
@click.option('--debug', is_flag=True)

def cmdline(path_cgc_original, path_cohorts, dict_mapping_cgc, dict_mapping_cgc_intogen, path_output,debug=False):
    bglogs.configure(debug=debug)
    cgc_dataset, intogen_ctypes, mapping, mapping_cgc_intogen=load_data(path_cgc_original,path_cohorts, dict_mapping_cgc, dict_mapping_cgc_intogen)
    cgc_dataset["acronym_cgc"] = map_cgc_into_acronym(cgc_dataset,mapping)
    intogen_ttypes = []

    for index,row in cgc_dataset.iterrows():
        intogen_ttypes.append([map_acronyms_into_intogen(row["acronym_cgc"],mapping_cgc_intogen)])
    cgc_dataset["cancer_type"] = ",".join(intogen_ttypes)
    # Save the output
    os.makedirs(path_output, exist_ok=True)
    output_file = os.path.join(path_output, 'cancer_gene_census_parsed.tsv')
    cgc_dataset.to_csv(output_file,sep="\t",index=False)


if __name__ == "__main__":
    cmdline()







