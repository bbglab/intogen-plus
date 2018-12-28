## IntOGen 2017 ##

# INSTALL

1. Install miniconda or conda python distributions
2. Create two environments:
    # conda create -n intogen2017 --file conda_env_intogen2017.txt
3. Edit configuration files at config/*.conf

# USAGE

nextflow run intogen.nf -profile local --input ./input --output ./output

