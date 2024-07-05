process FilterMNVS {
	tag "MNVs filter"
	publishDir "${STEPS_FOLDER}/boostDM", mode: "copy"
	label "core"

	input:
		path(input) from FILT_MNVS_INPUTS.collect()

	output:
		path("mnvs.tsv.gz") into MNVS_FILTER
	
	script:
		"""
		
		parse-mnvs ${input}

		"""
}