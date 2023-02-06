#!/usr/bin/env bash
set -e

FOLDER=$1
IMAGE=$2

tmpdir=`mktemp -d`

## if you have root permission you can uncomment the following:

# (cd ${FOLDER} && sudo singularity build ${tmpdir}/$(basename ${IMAGE}) Singularity)
# mv ${tmpdir}/$(basename ${IMAGE}) ${IMAGE}
# sudo chown ${USER}: ${IMAGE}

## if you DO NOT have root permission then you could use the --fakeroot flag as long as your user is specified in the /etc/subuid and /etc/subgid files. See https://docs.sylabs.io/guides/3.6/admin-guide/user_namespace.html#user-namespaces-fakeroot for more info. 

(cd ${FOLDER} && singularity build --fakeroot ${tmpdir}/$(basename ${IMAGE}) Singularity)
mv ${tmpdir}/$(basename ${IMAGE}) ${IMAGE}
chown ${USER}: ${IMAGE}
