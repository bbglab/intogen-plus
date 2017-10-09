#!/usr/bin/env nextflow

INPUT = file(params.input)
OUTPUT = file(params.output)


process PreprocessFromInput {
    tag { 'Reading variants' }
    publishDir OUTPUT, mode: 'copy'
    afterScript "cp .command.log $OUTPUT/preprocess_from_input.log"

    output:
        file "vep/*.in.gz" into IN_VEP mode flatten
        file "oncodrivefml/*.in.gz" into IN_ONCODRIVEFML mode flatten

    """
    python $baseDir/intogen4.py preprocess --cores $task.cpus -i $INPUT -o . vep oncodrivefml
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
        file "hotmapssignature/*.in.gz" into IN_HOTMAPS mode flatten
        file "oncodriveclust/*.in.gz" into IN_ONCODRIVECLUST mode flatten
        file "mutsigcv/*.in.gz" into IN_MUTSIGCV mode flatten

    """
    python $baseDir/intogen4.py read -i $task_file -o . oncodriveomega hotmapssignature oncodriveclust mutsigcv
    """
}

process OncodriveClust {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_ONCODRIVECLUST

    output:
        val task_file into OUT_ONCODRIVECLUST

    """
    if [ ! -f "${outputFile(OUTPUT, 'oncodriveclust', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . oncodriveclust $task_file
    fi
    """
}

process OncodriveFML {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_ONCODRIVEFML

    output:
        val task_file into OUT_ONCODRIVEFML

    """
    if [ ! -f "${outputFile(OUTPUT, 'oncodrivefml', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . oncodrivefml $task_file
    fi
    """
}

process HotmapsSignature {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_HOTMAPS

    output:
        val task_file into OUT_HOTMAPS

    """
    if [ ! -f "${outputFile(OUTPUT, 'hotmapssignature', task_file)}" ]
    then
        export PROCESS_CPUS=$task.cpus
        python $baseDir/intogen4.py run -o . hotmapssignature $task_file
    fi
    """
}


process OncodriveOmega {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_ONCODRIVEOMEGA

    output:
        val task_file into OUT_ONCODRIVEOMEGA

    """
    if [ ! -f "${outputFile(OUTPUT, 'oncodriveomega', task_file)}" ]
    then
        export PROCESS_CPUS=$task.cpus
        python $baseDir/intogen4.py run -o . oncodriveomega $task_file
    fi
    """
}

process MutsigCV {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_MUTSIGCV

    output:
        val task_file into OUT_MUTSIGCV

    """
    if [ ! -f "${outputFile(OUTPUT, 'mutsigcv', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . mutsigcv $task_file
    fi
    """
}


IN_SCHULZE = OUT_ONCODRIVECLUST
                    .phase(OUT_ONCODRIVEFML){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_HOTMAPS){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_ONCODRIVEOMEGA){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_MUTSIGCV){ it -> it.fileName }
                    .map{ it[0] }

process Schulze {
    tag { task_file.fileName }

    input:
        val task_file from IN_SCHULZE

    """
    if [ ! -f "${outputFile(OUTPUT, 'schulze', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o $OUTPUT schulze $task_file[0]
    fi
    """
}



def outputFile(output_folder, process_folder, task_file) {
    return output_folder.toString() + '/' + process_folder + '/' + task_file.fileName.toString().replace('.in.gz', '.out.gz')
}