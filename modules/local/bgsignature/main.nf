REGIONS_PREFIX = ['WXS': 'cds', 'WGS': 'wg']

process ComputeProfile {
	tag "ComputeProfile ${cohort}"
	label "bgsignature"
	publishDir "${STEPS_FOLDER}/signature", mode: "copy"

	input:
		tuple val(cohort), path(input), val(platform) from VARIANTS_SIG.join(PLATFORMS2)

	output:
		tuple val(cohort), path(output) into SIGNATURES

	script:
		prefix = REGIONS_PREFIX[platform]
		output = "${cohort}.sig.json"
		if (prefix)
			"""
			bgsignature normalize -m ${input} \
				-r ${params.datasets}/regions/${prefix}.regions.gz \
				--normalize ${params.datasets}/signature/${prefix}.counts.gz \
				-s 3 -g hg38 --collapse \
				--cores ${task.cpus} \
				-o ${output}
			"""
		else
			error "Invalid prefix. Check platform: $platform"

}
