#!/usr/bin/env nextflow

INPUT = file(params.input)
OUTPUT = file(params.output)


process PreprocessFromInput {
    tag { 'Reading variants' }
    publishDir OUTPUT, mode: 'copy'
    afterScript "cp .command.log $OUTPUT/preprocess_from_input.log"

    output:
        file "vep/*.in.gz" into IN_VEP mode flatten
        file "filters/*.json" into FILTERS_VARIANTS

    """
    python $baseDir/intogen4.py preprocess --cores $task.cpus -i $INPUT -o . vep
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
    if [ ! -f "${outputFile(OUTPUT, 'vep', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . vep $task_file
    else
        mkdir -p ./vep && cp ${outputFile(OUTPUT, 'vep', task_file)} ./vep/
    fi
    """
}

process PreprocessFromVep {
    tag { task_file.fileName }

    publishDir OUTPUT, mode: 'copy'
    afterScript "cp .command.log $OUTPUT/preprocess_from_vep.log"

    input:
        val task_file from OUT_VEP

    output:
        file "oncodriveomega/*.in.gz" into IN_ONCODRIVEOMEGA mode flatten
        file "mutsigcv/*.in.gz" into IN_MUTSIGCV mode flatten
        file "filters/vep/*.json" into FILTERS_VEP

    """
    python $baseDir/intogen4.py read -i $task_file -o . oncodriveomega mutsigcv
    """
}

process OncodriveOmega {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_ONCODRIVEOMEGA

    output:
        file "oncodriveomega/*.out.gz" into OUT_ONCODRIVEOMEGA mode flatten

    """
    if [ ! -f "${outputFile(OUTPUT, 'oncodriveomega', task_file)}" ]
    then
        export PROCESS_CPUS=$task.cpus
        python $baseDir/intogen4.py run -o . oncodriveomega $task_file
    else
        mkdir -p ./oncodriveomega && cp ${outputFile(OUTPUT, 'oncodriveomega', task_file)} ./oncodriveomega/
    fi
    """
}

process MutsigCV {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_MUTSIGCV

    output:
        file "mutsigcv/*.out.gz" into OUT_MUTSIGCV mode flatten

    """
    if [ ! -f "${outputFile(OUTPUT, 'mutsigcv', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . mutsigcv $task_file
    else
        mkdir -p ./mutsigcv && cp ${outputFile(OUTPUT, 'mutsigcv', task_file)} ./mutsigcv/
    fi
    """
}

def outputFile(output_folder, process_folder, task_file) {
    return output_folder.toString() + '/' + process_folder + '/' + task_file.fileName.toString().replace('.in.gz', '.out.gz')
}