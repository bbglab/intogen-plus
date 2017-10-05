#!/usr/bin/env nextflow

INPUT = file(params.input)
OUTPUT = file(params.output)


process PreprocessFromInput {
    tag { 'Reading variants' }
    publishDir OUTPUT

    output:
        file "vep/*.in.gz" into IN_VEP
        file "oncodrivefml/*.in.gz" into IN_ONCODRIVEFML

    """
    python $baseDir/intogen4.py preprocess -i $INPUT -o . vep oncodrivefml
    """

}

process Vep {
    tag { task_file.fileName }
    publishDir OUTPUT

    input:
        val task_file from IN_VEP

    output:
        file "vep/*.out.gz" into OUT_VEP
        file "vep/logs/*.log" into LOG_VEP

    """
    python $baseDir/intogen4.py run -o . vep $task_file
    """
}

process PreprocessFromVep {
    tag { task_file.fileName }
    publishDir OUTPUT

    input:
        val task_file from OUT_VEP

    output:
        file "oncodriveomega/*.in.gz" into IN_ONCODRIVEOMEGA
        file "hotmaps/*.in.gz" into IN_HOTMAPS

    """
    python $baseDir/intogen4.py read -i $task_file -o . oncodriveomega hotmaps
    """
}

process OncodriveFML {
    tag { task_file.fileName }
    publishDir OUTPUT

    input:
        val task_file from IN_ONCODRIVEFML

    output:
        file "oncodrivefml/*.out.gz" into OUT_ONCODRIVEFML
        file "oncodrivefml/logs/*.log" into LOG_ONCODRIVEFML

    """
    python $baseDir/intogen4.py run -o . oncodrivefml $task_file
    """
}

process OncodriveOmega {
    tag { task_file.fileName }
    publishDir OUTPUT

    input:
        val task_file from IN_ONCODRIVEOMEGA

    output:
        file "oncodriveomega/*.out.gz" into OUT_ONCODRIVEOMEGA
        file "oncodriveomega/logs/*.log" into LOG_ONCODRIVEOMEGA

    """
    export PROCESS_CPUS=$task.cpus
    python $baseDir/intogen4.py run -o . oncodriveomega $task_file
    """
}

process HotmapsSignature {
    tag { task_file.fileName }
    publishDir OUTPUT

    input:
        val task_file from IN_HOTMAPS

    output:
        file "hotmapssignature/*.out.gz" into OUT_HOTMAPS
        file "hotmapssignature/logs/*.log" into LOG_HOTMAPS

    """
    export PROCESS_CPUS=$task.cpus
    python $baseDir/intogen4.py run -o . hotmaps $task_file
    """
}