#!/usr/bin/env nextflow

INPUT = file(params.input)
OUTPUT = file(params.output)


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
    python $baseDir/intogen4.py preprocess --cores $task.cpus -i $INPUT -o . vep oncodrivefml dndscv oncodriveclustl smregions
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

    input:
        val task_file from OUT_VEP

    output:
        file "hotmapssignature/*.in.gz" into IN_HOTMAPS mode flatten
        file "edriver/*.in.gz" into IN_EDRIVER mode flatten
        file "cbase/*.in.gz" into IN_CBASE mode flatten
        file "filters/vep/*.json" into FILTERS_VEP
        file "deconstructsig/*.in.gz" into IN_DECONSTRUCTSIG mode flatten

    """
    python $baseDir/intogen4.py read -i $task_file -o . hotmapssignature edriver cbase deconstructsig
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
    if [ ! -f "${outputFile(OUTPUT, 'deconstructsig', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . deconstructsig $task_file
    else
        mkdir -p ./deconstructsig && cp -r ${outputFile(OUTPUT, 'deconstructsig', task_file)} ./deconstructsig/
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
        mkdir -p ./dndscv && cp ${outputPattern(OUTPUT, 'dndscv', task_file)} ./dndscv/
    fi
    """
}

OUT_DNDSCV
.filter({x -> !(x.fileName.toString().contains("genemuts") || x.fileName.toString().contains("annotmuts"))})
.into { OUT_DNDSCV_FOR_COMBINATION; OUT_DNDSCV_FOR_MUTRATE }

process MutRate {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'move'

    input:
        val task_file from OUT_DNDSCV_FOR_MUTRATE

    output:
        file "mutrate/*" into OUT_MUTRATE mode flatten

    """
    if [ ! -f "" ]
    then
        export INTOGEN_CPUS=$task.cpus
        python $baseDir/intogen4.py run -o . mutrate $task_file
    else
        mkdir -p ./mutrate && cp -r ${outputFile(OUTPUT, 'mutrate', task_file)} ./mutrate/
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
        export INTOGEN_CPUS=$task.cpus
        python $baseDir/intogen4.py run -o . oncodrivefml $task_file
    else
        mkdir -p ./oncodrivefml && cp ${outputFile(OUTPUT, 'oncodrivefml', task_file)} ./oncodrivefml/
    fi
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
    if [ ! -f "${outputFile(OUTPUT, 'smregions', task_file)}" ]
    then
        export INTOGEN_CPUS=$task.cpus
        python $baseDir/intogen4.py run -o . smregions $task_file
    else
        mkdir -p ./smregions && cp ${outputFile(OUTPUT, 'smregions', task_file)} ./smregions/
    fi
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
    if [ ! -f "${outputFile(OUTPUT, 'oncodriveclustl', task_file)}" ]
    then
        export INTOGEN_CPUS=$task.cpus
        export VEP_OUTPUT=$OUTPUT
        python $baseDir/intogen4.py run -o . oncodriveclustl $task_file
    else
        mkdir -p ./oncodriveclustl && cp ${outputFile(OUTPUT, 'oncodriveclustl', task_file)} ./oncodriveclustl/
    fi
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
    if [ ! -f "${outputFile(OUTPUT, 'hotmapssignature', task_file)}" ]
    then
        export INTOGEN_CPUS=$task.cpus
        python $baseDir/intogen4.py run -o . hotmapssignature $task_file
    else
        mkdir -p ./hotmapssignature && cp ${outputFile(OUTPUT, 'hotmapssignature', task_file)} ./hotmapssignature/
    fi
    """
}

process EDriver {
    tag { task_file.fileName }
    publishDir OUTPUT, mode: 'copy'

    input:
        val task_file from IN_EDRIVER

    output:
        file "edriver/*.out.gz" into OUT_EDRIVER mode flatten

    """
    if [ ! -f "${outputFile(OUTPUT, 'edriver', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . edriver $task_file
    else
        mkdir -p ./edriver && cp ${outputFile(OUTPUT, 'edriver', task_file)} ./edriver/
    fi
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
    if [ ! -f "${outputFile(OUTPUT, 'cbase', task_file)}" ]
    then
        python $baseDir/intogen4.py run -o . cbase $task_file
    else
        mkdir -p ./cbase && cp ${outputFile(OUTPUT, 'cbase', task_file)} ./cbase/
    fi
    """
}


IN_COMBINATION = OUT_ONCODRIVEFML
                    .phase(OUT_ONCODRIVECLUSTL){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_HOTMAPS){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_DNDSCV_FOR_COMBINATION){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_EDRIVER){ it -> it.fileName }
                    .map{ it[0] }
                    .phase(OUT_CBASE){ it -> it.fileName }
                    .map{ it[0] }

process Combination {
    tag { task_file.fileName }

    input:
        val task_file from IN_COMBINATION

    """
    if [ ! -f "${outputCombination(OUTPUT, task_file)}" ]
    then
        python $baseDir/intogen4.py run -o $OUTPUT combination $task_file
    fi
    """
}


def outputCombination(output_folder, task_file) {
    return output_folder.toString() + '/combination/' + task_file.fileName.toString().replace('.out.gz', '.05.out.gz')
}

def outputFile(output_folder, process_folder, task_file) {
    return output_folder.toString() + '/' + process_folder + '/' + task_file.fileName.toString().replace('.in.gz', '.out.gz')
}

def outputPattern(output_folder, process_folder, task_file) {
    return output_folder.toString() + '/' + process_folder + '/' + task_file.fileName.toString().replace('.in.gz', '*.out.gz')
}
