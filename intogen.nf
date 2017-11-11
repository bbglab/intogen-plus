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
        file "dndscv/*.in.gz" into IN_DNDSCV mode flatten
        file "filters/*.json" into FILTERS_VARIANTS

    """
    python $baseDir/intogen4.py preprocess --cores $task.cpus -i $INPUT -o . vep oncodrivefml dndscv
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
        file "hotmapssignature/*.in.gz" into IN_HOTMAPS mode flatten
        file "oncodriveclust/*.in.gz" into IN_ONCODRIVECLUST mode flatten
        file "filters/vep/*.json" into FILTERS_VEP

    """
    python $baseDir/intogen4.py read -i $task_file -o . hotmapssignature oncodriveclust
    """
}

process OncodriveClust {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_ONCODRIVECLUST

    output:
        file "oncodriveclust/*.out.gz" into OUT_ONCODRIVECLUST mode flatten

    """
    if [ ! -f "${outputFile(OUTPUT, 'oncodriveclust', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . oncodriveclust $task_file
    else
        mkdir -p ./oncodriveclust && cp ${outputFile(OUTPUT, 'oncodriveclust', task_file)} ./oncodriveclust/
    fi
    """
}

process DndsCV {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_DNDSCV

    output:
        file "dndscv/*.out.gz" into OUT_DNDSCV mode flatten

    """
    if [ ! -f "${outputFile(OUTPUT, 'dndscv', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . dndscv $task_file
    else
        mkdir -p ./dndscv && cp ${outputFile(OUTPUT, 'dndscv', task_file)} ./dndscv/
    fi
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
    if [ ! -f "${outputFile(OUTPUT, 'oncodrivefml', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . oncodrivefml $task_file
    else
        mkdir -p ./oncodrivefml && cp ${outputFile(OUTPUT, 'oncodrivefml', task_file)} ./oncodrivefml/
    fi
    """
}

process HotmapsSignature {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    maxForks 40

    input:
        val task_file from IN_HOTMAPS

    output:
        file "hotmapssignature/*.out.gz" into OUT_HOTMAPS mode flatten

    """
    if [ ! -f "${outputFile(OUTPUT, 'hotmapssignature', task_file)}" ]
    then
        export PROCESS_CPUS=$task.cpus
        python $baseDir/intogen4.py run -o . hotmapssignature $task_file
    else
        mkdir -p ./hotmapssignature && cp ${outputFile(OUTPUT, 'hotmapssignature', task_file)} ./hotmapssignature/
    fi
    """
}


IN_SCHULZE = OUT_ONCODRIVECLUST
                    .phase(OUT_ONCODRIVEFML){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_HOTMAPS){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_DNDSCV){ it -> it.fileName }
                    .map{ it[0] }

process Schulze {
    tag { task_file.fileName }

    input:
        val task_file from IN_SCHULZE

    """
    if [ ! -f "${outputFile(OUTPUT, 'schulze', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o $OUTPUT schulze $task_file
    fi
    """
}



def outputFile(output_folder, process_folder, task_file) {
    return output_folder.toString() + '/' + process_folder + '/' + task_file.fileName.toString().replace('.in.gz', '.out.gz')
}