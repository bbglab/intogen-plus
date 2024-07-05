process CBaSE {
	tag "CBaSE ${cohort}"
	publishDir "${STEPS_FOLDER}/cbase", mode: "copy"

    input:
        tuple val(cohort), path(input) from VARIANTS_CBASE

    output:
        tuple val(cohort), path(output) into OUT_CBASE

	script:
		output = "${cohort}.cbase.tsv.gz"
		"""
		mkdir -p Output/

		python /cbase/cbase.py ${input} ${params.datasets}/cbase 0 output
		tail -n+2 Output/q_values_output.txt | gzip > ${output}
		"""
}