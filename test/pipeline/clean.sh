#!/bin/bash

# Clean all NextFlow runs
rm -rf .nextflow
rm trace.txt*
rm timeline.html*
rm .nextflow.log*
ls work/ | grep -v conda | awk '{system("echo "$0" && rm -rf work/"$0)}'