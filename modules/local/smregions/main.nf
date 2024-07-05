process SMREGIONS {
	tag "SMRegions ${cohort}"
	publishDir "${STEPS_FOLDER}/smregions", mode: "copy"

    input:
        tuple val(cohort), path(input), path(signature)  from VARIANTS_SMREGIONS.join(SIGNATURES3)
        path regions from REGIONS

    output:
        tuple val(cohort), path(output) into OUT_SMREGIONS

	script:
		output = "${cohort}.smregions.tsv.gz"
		seedOpt = (params.seed == null)? '': "--seed ${params.seed}"
		debugOpt =  (params.debug)? '--debug': ''
		"""
		smregions -m ${input} -e ${regions} \
			-r ${params.datasets}/smregions/regions_pfam.tsv \
			-s ${signature} --cores ${task.cpus} \
			-c /smregions/smregions.conf \
			-o ${output} ${seedOpt} ${debugOpt}

		cat <<-END_VERSIONS > versions.yml
		"${task.process}":
			openvariant: \$(echo \$(openvar --version 2>&1) | sed 's/^.*openvar, version //; s/ *\$//' ))
		END_VERSIONS
		"""
}