## IntOGen ##

# INSTALL

1. Install miniconda or conda python distributions
2. Install singularity
3. Install nextflow
4. Edit configuration files at config/*.conf

# USAGE
cd datasets/test
nextflow run ../../intogen.nf -profile local --input ./input --output ./output

