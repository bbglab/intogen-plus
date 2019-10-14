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
<Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "chromosome_name" value = "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y"/>
        <Filter name = "biotype" value = "protein_coding"/>
		<Attribute name = "ensembl_gene_id" />
		<Attribute name = "external_gene_name" />
		<Attribute name = "ensembl_peptide_id" />
		<Attribute name = "chromosome_name" />
		<Attribute name = "genomic_coding_start" />
		<Attribute name = "genomic_coding_end" />
		<Attribute name = "cds_start" />
		<Attribute name = "cds_end" />
		<Attribute name = "cds_length" />
		<Attribute name = "strand" />
        <Attribute name = "ensembl_transcript_id" />
	</Dataset>
</Query>
'''
BIOMART_QUERY_ENCODED=$(python3 -c "from urllib.parse import quote_plus; print(quote_plus('''${BIOMART_QUERY}'''.replace('\n', '')))")

CDS_BIOMART_FILE="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/cds_biomart.tsv"
CDS_REGIONS_FILE="../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/cds.regions.gz"

mkdir -p ../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared
curl -s ${BIOMART_URL}?query=${BIOMART_QUERY_ENCODED} | grep -f <(cut -f2 ../../datasets/${INTOGEN_GENOME}_${INTOGEN_VEP}_${INTOGEN_RELEASE}/shared/ensembl_canonical_transcripts.tsv) | awk -F'\t' '($5!=""){print($0)}' > ${CDS_BIOMART_FILE}

cat ${CDS_BIOMART_FILE} | awk -F'\t' '($5!=""){gsub("-1", "-", $10); gsub("1", "+", $10); print($4"\t"$5"\t"$6"\t"$10"\t"$1"\t"$1"\t"$2)}' > "tmpfile"
echo -e "CHROMOSOME\tSTART\tEND\tSTRAND\tELEMENT\tSEGMENT\tSYMBOL" > "tmpheader"
cat "tmpheader" "tmpfile" | gzip > ${CDS_REGIONS_FILE}
rm "tmpheader" "tmpfile"
