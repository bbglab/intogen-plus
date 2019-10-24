Create the environment
----------------------

conda create -n intogen-datasets python=3.6
conda activate intogen-datasets
conda install mysql 
conda install bgsignature
conda install pandas
conda install tabix
conda install -c conda-forge -c bioconda ensembl-vep=92.4