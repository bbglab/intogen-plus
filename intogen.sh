#!/usr/bin/env bash
set -e


function run() {
	nextflow run intogen.nf "$@"
}

function test() {
	nextflow run test.nf --input ${PWD}/test/pipeline/input/cbioportal_prad_broad --output test/output --debug true -profile local -resume
}


function clean() {
	rm -r .nextflow*
	rm -r work
	rm trace.txt*
	rm timeline.html*
}

function usage()
{
    echo "-h | --help"
    echo ""
    echo "clean|run|test"
}


while (( "$#" ))
do
  case "$1" in
        -h | --help)
            usage
            exit
            ;;
        run)
            run "${@:1}"
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
        *)
            echo "ERROR: unknown parameter \"$1\""
            usage
            exit 1
            ;;
	esac
    shift
done
