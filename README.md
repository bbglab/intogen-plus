## IntOGen ##

# Install requirements

1. Install miniconda or conda python distributions
2. Install singularity
3. Install nextflow  (You can use `conda install nextflow`)

# Download and build datasets
All the datasets will be build at `datasets` folder. If you want to use a different folder create a symbolic link before running this script.

```bash
cd datasets_build
./build_all_datasets.sh hg38 vep92 develop
```

# Build singularity images
All the singularity images will be build at `containers` folder. If you want to use a different folder create a symbolic link before running this script.

```bash
cd containers_build
./build_all_containers.sh hg38 vep92 develop
```

# Test the pipeline
```bash
cd test/pipeline
conda activate
./run_hg38.sh
```
