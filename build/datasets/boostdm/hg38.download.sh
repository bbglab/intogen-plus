#!/bin/bash

set -xe

# define the paths
DEST=$1 

path_output=${DEST}/hg38.phyloP100way.bw

rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw  ${path_output}
