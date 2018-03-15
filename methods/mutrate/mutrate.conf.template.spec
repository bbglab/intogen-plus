[base]
# basic annotations
site_counts_file = string(default='/home/fmuinos/projects/intogen-plus/methods/mutrate/site_counts.pickle.gz')
triplets = string(default='/home/fmuinos/projects/masters/RB/regions/tnt_counts_dict.pickle.gz')
hugo_ensembl = string(default='/home/fmuinos/projects/intogen-plus/methods/mutrate/HGNC_ENSEMBL_dict.pickle')
ensembl_hugo = string(default='/home/fmuinos/projects/intogen-plus/methods/mutrate/ENSEMBL_HGNC_dict_70.pickle')

[sigfit]
# signature fitting data
cosmic = string(default='/home/fmuinos/projects/masters/RB/cosmic_signatures/cosmic_signatures.pickle.gz')
cosmic_exome = string(default='/home/fmuinos/projects/masters/RB/cosmic_signatures/cosmic_exome.tsv')