## IntOGen ##

# Install requirements

1. Install miniconda or conda python distributions
2. Install singularity
3. Install nextflow  (You can use `conda install nextflow`)

# Download or build datasets

All the required datasets will be stored at `datasets` folder,  if you want to use a different folder you must create a symbolic link before running this script. 

You can download all the required datasets using this command:

```bash
./download_datasets.sh hg38 vep92 stable
```

Or alternatively you can build all the datasets from the original sources. Note that this process can take a very long time and it might be broken if the original sources had changed.

```bash
cd datasets_build
./build_datasets.sh hg38 vep92 develop
```

# Download or build singularity images

All the images will be stored at `containers` folder, if you want to use a different folder you must create a symbolic link before running this script. 

You can download all the singularity images using this command:

```bash
./download_containers.sh hg38 vep92 stable
```

Or alternatively you can build all the singularity images yourself. 

```bash
cd containers_build
./build_containers.sh hg38 vep92 develop
```

# Test the pipeline
```bash
cd test/pipeline
conda activate
./run_hg38.sh
```
