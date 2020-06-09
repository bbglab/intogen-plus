# IntOGen #

> :warning: Please note that IntOGen needs a lot of resources. We strongly suggest to run it in a cluster environment!

### Install requirements

1. Install [miniconda or conda python distributions](https://docs.conda.io/projects/conda/en/latest/index.html)
2. Install [singularity](https://sylabs.io/singularity/) (the pipeline has been tested 
with version 2.x. You can use `conda install -c bioconda singularity`)
3. Install nextflow  (You can use `conda install nextflow`)
4. Clone this repository

### Download and build prerequisites

The IntOGen pipeline requires a collection of datasets and
Singularity containers in order to run.

> :warning: The pipeline has been tested with **hg38** and **vep92**

All the required datasets will be stored into the `datasets` folder, 
if you want to use a different folder you must create a symbolic link 
before running this script. 

First, make sure you have all the required libraries listed in `REQUIREMENTS.txt`

All the images will be stored into the `containers` folder, if you want to use a 
different folder you must create a symbolic link before running this script. 

Then you can build all the datasets from the original sources. 
Note that this process can take a very long time and it might fail if the 
original sources had changed.



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


### Licensing

IntoGEn uses a variety of software tools and datasets that
are released under a variety of licenses.

If you are using it for research/academic purposes it should be
fine. For commercial usage, you need to revise the license of the
different projects and datasets used.
