#!/bin/bash

set -e

export INTOGEN_RELEASE="latest"
mkdir -p ../containers/${INTOGEN_RELEASE}

# CBase
echo "Building cbase.simg"
cd cbase
sudo singularity build ../../containers/${INTOGEN_RELEASE}/cbase.simg Singularity
cd ..



