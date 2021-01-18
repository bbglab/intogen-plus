#!/usr/bin/env bash
set -e

# This helper script is meant to be used by BBGLab member only

SRC_FOLDER="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"


function run() {
	nextflow run {SRC_FOLDER}/intogen.nf --input "$@" --output intogen_{today} --debug true -profile bbglab -resume --containers {SRC_FOLDER}/containers --datasets {SRC_FOLDER}/datasets --annotations ${PWD}/annotations.txt
}

function test() {
	# Must be executed within this directory
	nextflow run intogen.nf --input ${PWD}/test/pipeline/input/cbioportal_prad_broad --output test/output --debug true -profile bbglab -resume
}


function clean() {
	rm -r .nextflow*
	rm -r work
	rm trace.txt*
	rm timeline.html*
}

function prepare() {
	folder=${1:-.}  # use current folder if not present
	mkdir -p ${folder}
	cp ${SRC_FOLDER}/config/annotations.txt ${folder}
	cp ${SRC_FOLDER}/intogen.sh ${folder}
	sed -e "s@{SRC_FOLDER}@${SRC_FOLDER}@g" -i ${folder}/intogen.sh
	today=`date +%Y%m%d`
	sed -e "s@{today}@${today}@g" -i ${folder}/intogen.sh
}

function usage()
{
    echo "-h | --help"
    echo ""
    echo "clean|run|test|prepare"
}


while (( "$#" ))
do
  case "$1" in
        -h | --help)
            usage
            exit
            ;;
        run)
            run "${@:2}"
            exit
            ;;
		test)
            test
            exit
            ;;
		clean)
			clean
			exit
			;;
		prepare)
            prepare "${@:2}"
            exit
            ;;
        *)
            echo "ERROR: unknown parameter \"$1\""
            usage
            exit 1
            ;;
	esac
    shift
done
