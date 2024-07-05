process dNdScv {
    tag "dNdScv ${cohort}"
    publishDir "${STEPS_FOLDER}/dndscv", mode: "copy"

    input:
        tuple val(cohort), path(input) from VARIANTS_DNDSCV

    output:
        tuple val(cohort), path("${cohort}.dndscv.tsv.gz") into OUT_DNDSCV
        tuple val(cohort), path("${cohort}.dndscv_annotmuts.tsv.gz") into OUT_DNDSCV_ANNOTMUTS
        tuple val(cohort), path("${cohort}.dndscv_genemuts.tsv.gz") into OUT_DNDSCV_GENEMUTS

	script:
		"""
		Rscript /dndscv/dndscv.R \
			${input} ${cohort}.dndscv.tsv.gz \
			${cohort}.dndscv_annotmuts.tsv.gz \
			${cohort}.dndscv_genemuts.tsv.gz
		"""
}