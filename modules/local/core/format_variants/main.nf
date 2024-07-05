process FormatSignature {
	tag "Prepare for signatures ${cohort}"
	label "core"
	publishDir "${STEPS_FOLDER}/signature", mode: "copy"

	input:
		tuple val(cohort), path(input) from VARIANTS1

	output:
		tuple val(cohort), path(output) into VARIANTS_SIG

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format ${type}
		"""

}