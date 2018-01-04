- Mutation data...

	genemuts.tsv
	annotmuts.tsv
	
	... were created using the following R script:
	
	#!/usr/bin/R
	
	library("seqinr")
	library("Biostrings")
	library("MASS")
	library("GenomicRanges")
	library("dndscv")
	data("dataset_simbreast", package="dndscv")
	dndsout = dndscv(mutations)
	write.table(dndsout$annotmuts, file="/home/fmuinos/projects/intogen-plus/methods/mutrate/annotmuts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
	write.table(dndsout$genemuts, file="/home/fmuinos/projects/intogen-plus/methods/mutrate/genemuts.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

- Site count data that contains site counts per gene, context and consequence...
	
	cds_site_count_dict_short.pickle.gz -- few consequence classes: 'syn', 'mis', 'trk';
	cds_site_count_dict.pickle.gz -- more comprehensive consequence classes

	... was copied from:

	/workspace/projects/oncodriveomega/site_counts/

