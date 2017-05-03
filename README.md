# Intogen 4


# Environment

Option 1

conda install variant-effect-predictor
vep_install.pl -a cf -s homo_sapiens -y GRCh37 -c /home/jordeu/workspace/intogen/intogen-plus/vepdata
vep_convert_cache.pl -species homo_sapiens -version 86_GRCh37 -d /home/jordeu/workspace/intogen/intogen-plus/vepdata


Option 2


conda install ensembl-vep
vep_install -a cf -s homo_sapiens -y GRCh38 -c /home/jordeu/workspace/intogen/intogen-plus/vepdata --CONVERT

