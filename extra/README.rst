
IntOGen pipeline extra content
==============================

In this document you can find scripts that we have used to do some
postprocessing of our data.
This is completely dependent of your own run, but here we provide instructions
in case you want to similar (or the same) things.

The Makefile will execute all this postprocessing of the data.


Requirements:

pytabix
numpy
pandas
bgoncotree

The following environment variables need to be defined:

- :envvar:`INTOGEN_INPUTS`: path to the info_datasets.csv files (folders that contain them)
- :envvar:`INTOGEN_RESULTS`: path to the output of the run of
  the pipeline

Optionally, this variable can also be defined:

- :envvar:`INTOGEN_DATASETS`: path to the datasets used in the
  pipeline. Defaults to :file:`../datasets`

