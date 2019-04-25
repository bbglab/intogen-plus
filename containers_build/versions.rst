IntOGen versions
================

Only versions that are specified during image build
are specified here.


Default
-------

This environment is the only one
that is used outside singularity


Base software:

    - Python (3.6)

Python libraries:

    - bgparsers
    - bgreference
    - click
    - intervaltree
    - numpy
    - pyliftover


CBaSE
-----

Base software: 

    - Python (2.7)

Python libraries:

    - mpmath
    - numpy
    - scipy
    - statsmodels

Custom scripts:

    - Python


Combination
-----------

Base software: 

    - Python (3.6)
    - bash

Python libraries:

    - click (7.0)
    - configobj (5.0)
    - cython (0.29)
    - numpy (1.15)
    - pandas (0.23)
    - scipy (1.2)
    - statsmodels (0.9)
    - tqdm (4.29)

Custom scripts:

    - bash
    - Python
    - configuration

DeconstructSig
--------------


Base software: 

    - Python (3)
    - R (3)

Python libraries:

    - click
    - pandas
    
R libraries:

    - deconstructSigs

.. note:: using BiocManager and BSgenome the strings
   for hg19 and hg38 genomes are installed as::

        BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

        BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")


Custom scripts:

    - Python
    - R
    - data


dNdScv
------

Base software: 

    - R (3)

R libraries:

    - dndscv (custom version)
         
Custom scripts:

    - R
    
HotMAPS
-------

Base software: 

    - Python (2.7)
    - bash
    - gzip

Python libraries:

    - bgreference
    - biopython
    - futures
    - pandas
    - pyliftover
    - tqdm
    
Custom scripts:

    - bash
    - Python


Mutrate
-------

Base software: 

    - Python (3)
    - bash
    - tar

Python libraries:

    - Click
    - numpy
    - pandas
    - tqdm
    
Custom scripts:

    - bash
    - Python
    
    
OncodriveCLUSTL
---------------


Base software: 

    - Python (3)

Python libraries:

    - oncodriveclustl
    
.. note:: currently we are using the version
   in the master branch of the repo
   
   

OncodriveFML
------------

Base software: 

    - Python (>=3.5)

Python libraries:

    - oncodrivefml
    
.. note:: we are using the *intogen* branch
   which is equivalent to v 2.2.0 but without
   generating the plots

    
Custom scripts:

    - configuration


Signature
---------

Base software: 

    - Python (3)
    - bash

Python libraries:

    - bgsignature
    
Custom scripts:

    - bash
   

SMRegions
---------

Base software: 

    - Python (3)

Python libraries:

    - smregions
    
.. note:: this package is provided
   with IntOGen source code

    
Custom scripts:

    - configuration

TransVar
--------


Base software: 

    - Python
    - samtools
    - tabix

Python libraries:

    - transvar (v2.4.3.20181231)



VEP
----

Base software: 

    - VEP (92.4)

