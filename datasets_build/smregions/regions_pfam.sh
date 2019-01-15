#!/bin/bash

set -e

if [ -z "${INTOGEN_RELEASE}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi

# Biomart release 92
BIOMART_URL="http://apr2018.archive.ensembl.org/biomart/martservice"

# Biomart Query
BIOMART_QUERY='''
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
		<Filter name = "biotype" value = "protein_coding"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "pfam_start" />
		<Attribute name = "pfam_end" />
		<Attribute name = "pfam" />
	</Dataset>
</Query>
'''
BIOMART_QUERY_ENCODED=$(python3 -c "from urllib.parse import quote_plus; print(quote_plus('''${BIOMART_QUERY}'''.replace('\n', '')))")

PFAM_BIOMART_FILE="../../datasets/${INTOGEN_RELEASE}/smregions/pfam_biomart.tsv"

mkdir -p ../../datasets/${INTOGEN_RELEASE}/smregions
curl -s ${BIOMART_URL}?query=${BIOMART_QUERY_ENCODED} | grep -f <(cut -f2 ../../datasets/${INTOGEN_RELEASE}/shared/ensembl_canonical_transcripts.tsv) | awk -F'\t' '($5!=""){print($0)}' > ${PFAM_BIOMART_FILE}



