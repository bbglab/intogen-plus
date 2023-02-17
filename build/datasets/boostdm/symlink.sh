#!/usr/bin/env bash
set -e 
DEST=$1
SOURCE=$2


echo "Creating symlink of ${SOURCE}"
ln -sfr ${SOURCE} ${DEST}

# python somescripts to produce  
# pfam_names.info.csv????
