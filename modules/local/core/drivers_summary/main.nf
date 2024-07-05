process DriverSummary {
	tag "Driver summary"
	publishDir "${OUTPUT}", mode: "copy"
	label "core"

    input:
        path (input) from DRIVERS.collect()
        path (input_vet) from VET.collect()
        path (mutations) from MUTATIONS_SUMMARY
        path (cohortsSummary) from COHORT_SUMMARY

    output:
		path("drivers.tsv") into DRIVERS_SUMMARY
		path("unique_drivers.tsv") into UNIQUE_DRIVERS
		path("unfiltered_drivers.tsv") into UNFILTER_DRIVERS

	script:
		"""
		drivers-summary \
			--mutations ${mutations} \
			--cohorts ${cohortsSummary} \
			${input} "${input_vet}"
		"""
}
