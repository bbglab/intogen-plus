# IntOGen #

> :warning: Please note that IntOGen needs a lot of resources. We strongly suggest to run it in a cluster environment!

### Install requirements

1. Install [miniconda or conda python distributions](https://docs.conda.io/projects/conda/en/latest/index.html)
2. Install [singularity](https://sylabs.io/singularity/) (the pipeline has been tested 
with version 2.x. You can use `conda install -c bioconda singularity`)
3. Install nextflow  (You can use `conda install nextflow`)

### Download and build datasets

> :warning: The pipeline has been tested with **hg38** and **vep92**

All the required datasets will be stored into the `datasets` folder, 
if you want to use a different folder you must create a symbolic link 
before running this script. 

First, make sure you have all the required libraries listed in `REQUIREMENTS.txt`

Then you can build all the datasets from the original sources. 
Note that this process can take a very long time and it might fail if the 
original sources had changed.

```bash
cd datasets_build
./build_datasets.sh hg38 vep92 stable
```

where:
```bash
INTOGEN_GENOME is the release of the genome (such as hg38)
INTOGEN_VEP is the release of VEP (such as vep92)
INTOGEN_RELEASE is the release of the intogen build (such as develop or stable)
```




### Download and build singularity images

> :warning: The pipeline has been tested with **hg38** and **vep92**

All the images will be stored into the `containers` folder, if you want to use a 
different folder you must create a symbolic link before running this script. 

You can build all the singularity images using this command:

```bash
cd containers_build
./build_containers.sh <INTOGEN_GENOME> <INTOGEN_VEP> <INTOGEN_RELEASE>
```

where:
```bash
INTOGEN_GENOME is the release of the genome (such as hg19, hg38)
INTOGEN_VEP is the release of VEP (such as vep88, vep92)
INTOGEN_RELEASE is the release of the intogen build (such as develop or stable)
```

### Run the pipeline

The IntOGen pipeline is built on top of [nextflow](https://www.nextflow.io/). 
In order to execute the pipeline, first activate the virtual environment 
you previously created with the required libraries. Then create a bash file 
with the following content but substituting the variables `INTOGEN_GENOME`, 
`INTOGEN_VEP`, `INTOGEN_RELEASE`, and `INTOGEN_HOME` with your settings:
```bash
#!/bin/bash

export INTOGEN_GENOME="hg38"
export INTOGEN_VEP="vep92"
export INTOGEN_RELEASE="stable"
export INTOGEN_HOME=`echo "$(pwd)/../../" | xargs realpath`

nextflow run ${INTOGEN_HOME}/intogen.nf -resume -profile local --input ./input --output ./output
```
If you want to execute the pipeline in a cluster environment, change the 
value of the nextflow argument `-profile` to `slurm`. Also, if you want 
to change the number of cores and memory assigned to each nextflow process,
you can edit the `local.conf` (if you set -profile local) or the `slurm.conf` 
(if you set -profile slurm) files in the `config` folder.

### Test the pipeline
The `test` folder include a dataset to test the pipeline: 
```bash
cd test/pipeline
conda activate <NAME_OF_YOUR_VIRTUALENVIRONMENT>
./run_hg38.sh
```