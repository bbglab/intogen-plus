process OPENVARIANT_CAT {
    tag "$meta.id"
    label 'openvariant'

    container "docker.io/federicabrando/openvariant:v0.7.6"

    input:
        tuple val(meta), path(input)
    output:
        tuple val(meta), path("*.parsed.tsv.gz")    , emit: parsed_output
        path ("versions.yml")                       , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cohort = input.baseName.split('\\.')[0]
    
    """
    openvar cat ${input} --header | gzip > ${cohort}.parsed.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        openvariant: \$(echo \$(openvar --version 2>&1) | sed 's/^.*openvar, version //; s/ *\$//' ))
    END_VERSIONS
    """
}