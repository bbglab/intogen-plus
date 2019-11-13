#!/usr/bin/env nextflow

INPUT = file(params.input)
OUTPUT = file(params.output)
INTOGEN_SCRIPT = "python $baseDir/intogen.py"


process PreprocessFromInput {
    tag { 'Reading variants' }
    publishDir OUTPUT, mode: 'copy'

    output:
        file "vep/*.in.gz" into IN_VEP mode flatten
        file "oncodrivefml/*.in.gz" into IN_ONCODRIVEFML mode flatten
        file "oncodriveclustl/*.in.gz" into IN_ONCODRIVECLUSTL mode flatten
        file "dndscv/*.in.gz" into IN_DNDSCV mode flatten
        file "signature/*.in.gz" into IN_SIGNATURES mode flatten
        file "filters/*.json" into FILTERS_VARIANTS

    """
    $INTOGEN_SCRIPT readvariants --cores $task.cpus -i $INPUT -o $OUTPUT vep oncodrivefml dndscv oncodriveclustl signature
    """
}

process Signature {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_SIGNATURES

    output:
        file "signature/*.json" into OUT_SIGNATURE mode flatten

    """
    $INTOGEN_SCRIPT calculatesignature --cores $task.cpus -o $OUTPUT signature $task_file
    """
}

// Duplicate this stream
OUT_SIGNATURE.into { OUT_SIGNATURE_01; OUT_SIGNATURE_02; OUT_SIGNATURE_03; OUT_SIGNATURE_04; }

process Vep {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_VEP

    output:
        file "vep/*.out.gz" into OUT_VEP mode flatten

    """
    $INTOGEN_SCRIPT run -o $OUTPUT vep $task_file
    """
}

// Duplicate this stream
OUT_VEP.into { OUT_VEP_01; OUT_VEP_02; }

process PreprocessFromVep {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from OUT_VEP_01

    output:
        file "hotmaps/*.in.gz" into IN_HOTMAPS mode flatten
        file "cbase/*.in.gz" into IN_CBASE mode flatten
        file "filters/vep/*.json" into FILTERS_VEP
        file "mutpanning/*.in.gz" into IN_MUTPANNING mode flatten
        file "deconstructsig/*.in.gz" into IN_DECONSTRUCTSIG mode flatten

    """
    $INTOGEN_SCRIPT readvep -i $task_file -o $OUTPUT hotmaps cbase mutpanning deconstructsig
    """
}

process PreprocessFromVepNonSynonymous {
    tag { task_file.fileName }

    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from OUT_VEP_02

    output:
        file "smregions/*.in.gz" into IN_SMREGIONS mode flatten

    """
    $INTOGEN_SCRIPT readvepnonsynonymous -i $task_file -o $OUTPUT smregions
    """
}

process DeconstructSig {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_DECONSTRUCTSIG

    output:
        file "deconstructsig/*.out.gz" into OUT_DECONSTRUCTSIG mode flatten

    """
    $INTOGEN_SCRIPT run -o $OUTPUT deconstructsig $task_file
    """
}

process DndsCV {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_DNDSCV

    output:
        file "dndscv/*.out.gz" into OUT_DNDSCV mode flatten
        file "dndscv/*.annotmuts.gz" into ANNOTMUTS_DNDSCV mode flatten
        file "dndscv/*.genemuts.gz" into GENEMUTS_DNDSCV mode flatten

    """
    $INTOGEN_SCRIPT run -o $OUTPUT dndscv $task_file
    """
}

// Duplicate this stream
OUT_DNDSCV.into { OUT_DNDSCV_01; OUT_DNDSCV_02; }

// Combination stream
// IN_MUTRATE = OUT_DNDSCV_01.phase(OUT_DECONSTRUCTSIG){ it -> it.fileName }.map{ it[0] }
IN_MUTRATE = OUT_DNDSCV_01.phase(OUT_DECONSTRUCTSIG){ it -> it.fileName.toString().tokenize('.').get(0) }.map{ it[0] }

process MutRate {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_MUTRATE

    output:
        file "mutrate/*.tar.gz" into OUT_MUTRATE mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT mutrate $task_file
    """
}

// Combination stream
IN_ONCODRIVEFML_PHASED = IN_ONCODRIVEFML.phase(OUT_SIGNATURE_01){ it -> it.fileName.toString().tokenize('.').get(0) }.map{ it[0] }

process OncodriveFML {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_ONCODRIVEFML_PHASED

    output:
        file "oncodrivefml/*.out.gz" into OUT_ONCODRIVEFML mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT oncodrivefml $task_file
    """
}

process MutPanning {

    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_MUTPANNING

    output:
        file "mutpanning/*.out.gz" into OUT_MUTPANNING mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT mutpanning $task_file
    """
}

// Combination stream
IN_SMREGIONS_PHASED = IN_SMREGIONS.phase(OUT_SIGNATURE_02){ it -> it.fileName.toString().tokenize('.').get(0) }.map{ it[0] }

process SMRegions {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_SMREGIONS_PHASED

    output:
        file "smregions/*.out.gz" into OUT_SMREGIONS mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT smregions $task_file
    """
}

// Combination stream
IN_ONCODRIVECLUSTL_PHASED = IN_ONCODRIVECLUSTL.phase(OUT_SIGNATURE_03){ it -> it.fileName.toString().tokenize('.').get(0) }.map{ it[0] }

process OncodriveClustl {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_ONCODRIVECLUSTL_PHASED

    output:
        file "oncodriveclustl/*.out.gz" into OUT_ONCODRIVECLUSTL mode flatten
        file "oncodriveclustl/*.clusters.gz" into CLUSTERS_ONCODRIVECLUSTL mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT oncodriveclustl $task_file
    """
}

// Combination stream
IN_HOTMAPS_PHASED = IN_HOTMAPS.phase(OUT_SIGNATURE_04){ it -> it.fileName.toString().tokenize('.').get(0) }.map{ it[0] }

process Hotmaps {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_HOTMAPS_PHASED

    output:
        file "hotmaps/*.out.gz" into OUT_HOTMAPS mode flatten
        file "hotmaps/*.clusters.gz" into CLUSTERS_HOTMAPS mode flatten
        file "hotmaps/*.tmp" optional true into TMP_OUTPUTS mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT hotmaps $task_file
    """
}

process CBase {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_CBASE

    output:
        file "cbase/*.out.gz" into OUT_CBASE mode flatten

    """
    $INTOGEN_SCRIPT run -o $OUTPUT cbase $task_file
    """
}

// Combination stream
IN_COMBINATION = OUT_ONCODRIVEFML
                    .phase(OUT_ONCODRIVECLUSTL){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_HOTMAPS){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_DNDSCV_02){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_SMREGIONS){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_MUTPANNING){ it -> it.fileName }
                    .map{ it[0] }
                    // .phase(OUT_CBASE){ it -> it.fileName }
                    // .map{ it[0] }


process Combination {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'move'

    input:
        val task_file from IN_COMBINATION

    """
    $INTOGEN_SCRIPT run -o $OUTPUT combination $task_file
    """
}

