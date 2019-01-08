## IntOGen ##

# INSTALL

1. Install miniconda or conda python distributions
2. Install singularity
3. Create intogen environment:
    # conda create -n intogen python=3.5 nextflow
4. Activate the environment
    # conda activate intogen
5. Update nextflow
    # nextflow self-update  
6. Install dependencies
    # conda install -c bbglab bgparsers click bgreference
7. Edit configuration files at config/*.conf

# USAGE
cd datasets/test
nextflow run ../../intogen.nf -profile local --input ./input --output ./output

