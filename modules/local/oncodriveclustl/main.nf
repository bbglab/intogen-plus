process OncodriveCLUSTL {
    tag "OncodriveCLUSTL ${cohort}"
    publishDir "${STEPS_FOLDER}/oncodriveclustl", mode: "copy"

    input:
        tuple val(cohort), path(input), path(signature), val(cancer) from VARIANTS_CLUSTL.join(SIGNATURES2).join(CANCERS1)
        path regions from REGIONS

    output:
        tuple val(cohort), path("${cohort}.elements_results.txt") into OUT_ONCODRIVECLUSTL
        tuple val(cohort), path("${cohort}.clusters_results.tsv") into OUT_ONCODRIVECLUSTL_CLUSTERS

	script:
		seedOpt = (params.seed == null)? '': "--seed ${params.seed}"
		debugOpt =  (params.debug)? '--log-level debug': ''
		if (['CM', 'SBCC', 'SSCC'].contains(cancer))
			"""
			oncodriveclustl -i ${input} -r ${regions} \
				-g hg38 -sim region_restricted -n 1000 -kmer 5 \
				-sigcalc region_normalized \
				--concatenate \
				-c ${task.cpus} \
				-o ${cohort} ${seedOpt} ${debugOpt}
			
			mv ${cohort}/elements_results.txt ${cohort}.elements_results.txt
			mv ${cohort}/clusters_results.tsv ${cohort}.clusters_results.tsv
			"""
		else
			"""
			oncodriveclustl -i ${input} -r ${regions} \
				-g hg38 -sim region_restricted -n 1000 -kmer 3 \
				-sig ${signature} --concatenate \
				-c ${task.cpus} \
				-o ${cohort} ${seedOpt} ${debugOpt}

			mv ${cohort}/elements_results.txt ${cohort}.elements_results.txt
			mv ${cohort}/clusters_results.tsv ${cohort}.clusters_results.tsv	
			"""
}