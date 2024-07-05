process DriverDiscovery {
	tag "Driver discovery ${cohort}"
	publishDir "${STEPS_FOLDER}/drivers", mode: "copy"
	label "core"

    input:
        tuple val(cohort), path(combination), path(deconstruct_in), path(sig_likelihood), path(smregions), path(clustl_clusters), path(hotmaps_clusters), path(dndscv), val(cancer) from OUT_COMBINATION.join(VARIANTS_DECONSTRUCTSIGS2).join(OUT_DECONSTRUCTSIGS_SIGLIKELIHOOD).join(OUT_SMREGIONS2).join(OUT_ONCODRIVECLUSTL_CLUSTERS).join(OUT_HOTMAPS_CLUSTERS).join(OUT_DNDSCV2).join(CANCERS3)

    output:
		path(output_drivers) into DRIVERS
		path(output_vet) into VET

	script:
		output_drivers = "${cohort}.drivers.tsv"
		output_vet = "${cohort}.vet.tsv"
		"""
		drivers-discovery --output_drivers ${output_drivers} \
			--output_vet ${output_vet} \
			--combination ${combination} \
			--mutations ${deconstruct_in} \
			--sig_likelihood ${sig_likelihood} \
			--smregions ${smregions} \
			--clustl_clusters ${clustl_clusters} \
			--hotmaps ${hotmaps_clusters} \
			--dndscv ${dndscv} \
			--ctype ${cancer} \
			--cohort ${cohort}
		"""
}