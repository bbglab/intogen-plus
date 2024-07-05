process ONCODRIVEFML {
    tag "$meta.id"
    publishDir "${STEPS_FOLDER}/oncodrivefml", mode: "copy"

	container 'docker.io/bbglab/oncodrivefml:latest'
    
	input:
        tuple val(meta), val(cohort)
		tuple val(meta), path(input)		 
		tuple val(meta), path(signature)  
        path regions from REGIONS
		
	output:
		tuple val(meta),
		path ("versions.yml")

	when:
    task.ext.when == null || task.ext.when
    
	output:
        tuple val(cohort), path("out/*.tsv.gz") into OUT_ONCODRIVEFML

	script:
		seedOpt = (params.seed == null)? '': "--seed ${params.seed}"
		debugOpt =  (params.debug)? '--debug': ''
		"""
		oncodrivefml -i ${input} -e ${regions} --signature ${signature} \
			-c /oncodrivefml/oncodrivefml_v2.conf  --cores ${task.cpus} \
			-o out ${seedOpt} ${debugOpt}
		
	cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrivefml: \$(echo \$(oncodrivefml --version | rev | cut -d ' ' -f1 | rev))
    END_VERSIONS
		"""
	
	stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        oncodrivefml: \$(echo \$(oncodrivefml --version | rev | cut -d ' ' -f1 | rev))
    END_VERSIONS
    """
}