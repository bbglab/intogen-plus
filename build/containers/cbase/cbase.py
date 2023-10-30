#!/usr/bin/env python

# ************************************************************
# * Cancer Bayesian Selection Estimation (CBaSE):			 *
# * Code accompanying Weghorn & Sunyaev, Nat. Genet. (2017). *
# *															 *
# * Author:		Donate Weghorn								 *
# *															 *
# * Copyright:	(C) 2017-2019 Donate Weghorn				 *
# *															 *
# * License:		Public Domain							 *
# *															 *
# * Version:		1.1										 *
# ************************************************************


# This version incorporates the following changes:
# - for sake of computation, this version skips the negative selection analysis;
# - this version changes the output format accordingly to display positive selection analysis;
# - add a selection clause for statistical model of choice based on the total mutation burden;
# - add an implementation of the Benjamini-Hochberg FDR method to compute q-values.

import sys
import subprocess


if len(sys.argv) < 5:
    print("Usage: python cbase.py <input_file> <auxiliary_path> <c_mode> <outname>")
    print("Arguments:")
    print("<input_file>: Somatic mutation data input file")
    print("<auxiliary_path>: Path to auxiliary input files folder")
    print("<c_mode>: 0 = trinucleotides, 1 = pentanucleotides")
    print("<outname>: Name trunk for labeling files")
    sys.exit(1)

#************************************** COMMAND LINE ARGS ***************************************

infile		= str(sys.argv[1])		#	somatic mutation data input file
filepath	= str(sys.argv[2])		#	path to auxiliary input files folder
c_mode		= int(sys.argv[3])		#	0 = trinucleotides, 1 = pentanucleotides
# mod_C		= int(sys.argv[4])		#	model choice: 0=all, 1=G, 2=IG, 3=EmixG, 4=EmixIG, 5=GmixG, 6=GmixIG
outname		= str(sys.argv[4])		#	name trunk for labeling files


subprocess.call(["python", "-W", "ignore", "Auxiliary/CBaSE_v1.1_parameters.py", "%s" %infile, "%s" %filepath, "%s" %c_mode, "%s" %outname])
subprocess.call(["python", "-W", "ignore", "Auxiliary/CBaSE_v1.1_qvalues.py", "%s" %outname])
