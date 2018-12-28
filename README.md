## IntOGen 2017 ##

# INSTALL

1. Install miniconda or conda python distributions
2. Install singularity
3. Create intogen environment:
    # conda create -n intogen python=3.5 nextflow
4. Activate the environment
    # conda activate intogen
5. Update nextflow
    # nextflow self-update    
6. Edit configuration files at config/*.conf

# USAGE
nextflow run intogen.nf -profile local --input ./input --output ./output

