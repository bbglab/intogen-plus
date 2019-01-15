#!/bin/bash

set -e

if [ -z "${INTOGEN_RELEASE}" ]
then
      echo "ERROR: Define the INTOGEN_RELEASE variable"
      exit -1
fi

#ENSEMBL_DATABASE="homo_sapiens_core_92_38"
ENSEMBL_DATABASE="homo_sapiens_core_75_37"

mkdir -p ../../datasets/${INTOGEN_RELEASE}/shared
SQL_QUERY="SELECT g.stable_id, t.stable_id, x.display_label FROM gene g JOIN transcript t ON (g.canonical_transcript_id = t.transcript_id) JOIN xref x ON (g.display_xref_id = x.xref_id AND g.biotype='protein_coding') LEFT JOIN external_db ed USING (external_db_id) WHERE ed.db_name = 'HGNC'"
mysql -u anonymous -h ensembldb.ensembl.org --column-names=FALSE -e "${SQL_QUERY};" ${ENSEMBL_DATABASE} > ../../datasets/${INTOGEN_RELEASE}/shared/ensembl_canonical_transcripts.tsv