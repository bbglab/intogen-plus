process ParseProfile {
	tag "Parsing profile ${cohort}"
	publishDir "${STEPS_FOLDER}/boostDM/mutrate", mode: "copy"
	label "core"

    input:
        tuple val(cohort), path(signature) from SIGNATURES5

    output:
		tuple val(cohort), path("*.mutrate.json") into OUT_MUTRATE

	script:
		"""
		parse-profile -i ${signature} -o ${cohort}
		"""

}