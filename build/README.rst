
Build
=====

This directory contains a Makefile to build
the datasets and containers (as
`Singularity <https://sylabs.io/docs/>`_
images).

Although you can find many ``*.mk`` files,
they are meant to be called form the Makefile in this directory
and not directly, as paths, or variables might not be properly defined.
However, you can get an idea of what is the build process for each dataset
or container.

Currently there are several variables that can be defined:

- ``ENSEMBL``: specifies the ensembl version
- ``CADD``: specifies the CADD version
- ``DATASETS``: path where to store the datasets
- ``CONTAINERS``: path where to store the containers
- ``CORES``: amount of cores to use by processes that can be parallelized

In addition, you can make use of other Make features, such as
parallel build:

.. code:: bash

	$ make -j 4

Moreover, there is a special target (``sudo-containers``) that
can be used to build the containers that require superuser privileges.


Requirements
************

Requires a working Internet connection
and the following software:

awk
cut
xz
curl
make
mysql
sqlite3
singularity (version 2.6.1)
tabix
python
	bgdata
	bgsignature
	bgreference
	bgvep
	click
	numpy
	pandas
	tqdm

This software (except singularity) can be installed with
`conda <https://docs.conda.io/en/latest/>`_.


.. important:: The MutPanning jar file needs to be downloaded
   from their website (http://www.cancer-genes.org/ )
   extracted and placed in the ``containers_build/mutpanning`` folder.
