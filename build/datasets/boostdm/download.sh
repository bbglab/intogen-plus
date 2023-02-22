FOLDER=$1

tmpdir=`mktemp -d`

wget -c https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam35.0/Pfam-A.clans.tsv.gz \
		-O ${tmpdir}/pfam_info.name.tsv.gz

gzip -d ${tmpdir}/pfam_info.name.tsv.gz
mv ${tmpdir}/pfam_info.name.tsv ${FOLDER}/pfam_info.name.tsv