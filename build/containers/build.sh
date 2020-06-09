#!/usr/bin/env bash
set -e

FOLDER=$1
IMAGE=$2

tmpdir=`mktemp -d`

(cd ${FOLDER} && sudo singularity build ${tmpdir}/$(basename ${IMAGE}) Singularity)
mv ${tmpdir}/$(basename ${IMAGE}) ${IMAGE}
sudo chown ${USER}: ${IMAGE}
