process deconstructSigs {
	tag "deconstructSigs ${cohort}"
	publishDir "${STEPS_FOLDER}/deconstructSigs", mode: "copy"

    input:
        tuple val(cohort), path(input) from VARIANTS_DECONSTRUCTSIGS1

    output:
        tuple val(cohort), path(output) into OUT_DECONSTRUCTSIGS
        tuple val(cohort), path("*.signature_likelihood") into OUT_DECONSTRUCTSIGS_SIGLIKELIHOOD

	script:
		output = "${cohort}.deconstructsigs.tsv.gz"
		likelihood = "${cohort}.signature_likelihood"
		"""
		python3 /deconstructsig/run_deconstruct.py \
			--input-file ${input} --weights ${output} \
			--build hg38
		python3 /deconstructsig/signature_assignment.py \
			--input-file ${output} \
			--output-file ${likelihood}
		"""
}