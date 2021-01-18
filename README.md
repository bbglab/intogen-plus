# IntOGen #

> :warning: Please note that IntOGen needs a lot of resources. We strongly suggest to run it in a cluster environment!

### Install requirements

1. Install [singularity](https://sylabs.io/singularity/) (the pipeline has been tested 
with version 2.x)
2. Install nextflow  (You can use `conda install nextflow`)
3. Clone this repository

### Download and build prerequisites

The IntOGen pipeline requires a collection of datasets and
Singularity containers in order to run.

> :warning: The pipeline has been tested with **hg38** and **vep92** and **vep101**

See the README.rst file in the ``build`` folder for further details.

Then you can build all the datasets from the original sources. 
Note that this process can take a very long time and it might fail if the 
original sources had changed.



### Run the pipeline

The IntOGen pipeline is built on top of [nextflow](https://www.nextflow.io/). 
In order to execute the pipeline, you only need to execute:

[comment]: <> (FIXME add example in test)

```bash
nextflow run intogen.nf -resume -profile local --input test/ --output ./output
```

For further details, please check our documentation: http://intogen.rtfd.io/

[comment]: <> (FIXME add example in test)

To avoid stopping the pipeline execution for one or a few incorrect
input, we have decided to ignore the errors of the steps by default.
However, we advise to review each of them carefully to 
understand the underlying reasons.


### Licensing

IntoGEn uses a variety of software tools and datasets that
are released under a variety of licenses.
In order to accommodate all of them, the pipeline itself is released
under GNU General Public License version 3.

If you are using it for research/academic purposes it should be
fine. For commercial usage, you need to revise the license of the
different software pieces and datasets used, as some of them
(e.g. CADD or CGC) are restricted for commercial usage.
