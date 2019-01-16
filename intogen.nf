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
        file "smregions/*.in.gz" into IN_SMREGIONS mode flatten
        file "filters/*.json" into FILTERS_VARIANTS

    """
    $INTOGEN_SCRIPT readvariants --cores $task.cpus -i $INPUT -o . vep oncodrivefml dndscv oncodriveclustl smregions
    """

}

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

process PreprocessFromVep {
    tag { task_file.fileName }

    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from OUT_VEP

    output:
        file "hotmapssignature/*.in.gz" into IN_HOTMAPS mode flatten
        file "cbase/*.in.gz" into IN_CBASE mode flatten
        file "filters/vep/*.json" into FILTERS_VEP
        file "deconstructsig/*.in.gz" into IN_DECONSTRUCTSIG mode flatten

    """
    $INTOGEN_SCRIPT readvep -i $task_file -o . hotmapssignature cbase deconstructsig
    """
}

process DeconstructSig {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_DECONSTRUCTSIG

    output:
        file "deconstructsig/*.out.gz" into OUT_DECONSTRUCTSIG

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

OUT_DNDSCV.into { OUT_DNDSCV_FOR_COMBINATION; OUT_DNDSCV_FOR_MUTRATE }

process MutRate {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'move'

    input:
        val task_file from OUT_DNDSCV_FOR_MUTRATE

    output:
        file "mutrate/*" into OUT_MUTRATE mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT mutrate $task_file
    """
}

process OncodriveFML {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_ONCODRIVEFML

    output:
        file "oncodrivefml/*.out.gz" into OUT_ONCODRIVEFML mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT oncodrivefml $task_file
    """
}

process SMRegions {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_SMREGIONS

    output:
        file "smregions/*.out.gz" into OUT_SMREGIONS mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT smregions $task_file
    """
}


process OncodriveClustl {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_ONCODRIVECLUSTL

    output:
        file "oncodriveclustl/*.out.gz" into OUT_ONCODRIVECLUSTL mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT oncodriveclustl $task_file
    """
}


process HotmapsSignature {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_HOTMAPS

    output:
        file "hotmapssignature/*.out.gz" into OUT_HOTMAPS mode flatten

    """
    $INTOGEN_SCRIPT run -c $task.cpus -o $OUTPUT hotmapssignature $task_file
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


IN_COMBINATION = OUT_ONCODRIVEFML
                    .phase(OUT_ONCODRIVECLUSTL){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_HOTMAPS){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_DNDSCV_FOR_COMBINATION){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_CBASE){ it -> it.fileName }
                    .map{ it[0] }

process Combination {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'move'

    input:
        val task_file from IN_COMBINATION

    output:
        file "combination/*.out.gz"

    """
    $INTOGEN_SCRIPT run -o $OUTPUT combination $task_file
    """
}