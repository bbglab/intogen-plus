
// TODO add in the docs the option to set here a list of values
// Set here a list of files or directories to use as ["/path/*", "/path2/file"]
INPUT = Channel.fromPath(params.input)
//INPUT = Channel.fromPath(["test/*", "test/pipeline"], type: 'any')
OUTPUT = file(params.output)
DEBUG_FOLDER = file(params.debugFolder)
// TODO add in the docs
ANNOTATIONS = Channel.value(params.annotations)
REGIONS = Channel.value("${params.datasets}/regions/cds.regions.gz")

/* TODO remove
process test {
	input:
		path p from INPUT


	script:
		"""
		echo $BGDATA_LOCAL
		"""
}*/


process parseInput {
	tag "Parse input ${input}"
	label "core"
	publishDir "${DEBUG_FOLDER}/inputs", mode: "symlink", enabled: params.debug

	input:
		path input from INPUT
		path annotations from ANNOTATIONS

	output:
	// TODO remove .parsed if possible
		// mode flatten?
		path("*.parsed.tsv.gz") into PARSED_INPUT

	script:
		// TODO only bginfo
		if ( input.toRealPath().toFile().isDirectory() || input.endsWith(".bginfo" ))
			// TODO explain that dataset is required
			"""
			bgvariants groupby --cores ${task.cpus} -s 'gzip > \${GROUP_KEY}.parsed.tsv.gz' --headers -g DATASET -a ${annotations} ${input}
			"""

		else
			// filename used as dataset name
			// TODO: create a .tsv.gz with the output and names as DATASET
			"""
			bgvariants cat -a ${annotations} ${input}
			"""
}

PARSED_INPUT.flatten().into{ PARSED_INPUT1; PARSED_INPUT2; PARSED_INPUT3 }


CUTOFFS = ['WXS': 1000, 'WGS': 10000]

process processVariants {
	tag "Process variants ${cohort}"
	label "core"
	errorStrategy 'ignore'  // if a cohort does not pass the filters, do not proceed with it
	publishDir "${DEBUG_FOLDER}/variants", mode: "symlink", enabled: params.debug

	input:
		path input from PARSED_INPUT1

	output:
		path output into VARIANTS
		tuple val(cohort), path("${output}.stats.json") into STATS_VARIANTS
		// FIXME create separate channels for platform and cancer type
		tuple val(cohort), val(platform) into PLATFORMS

	script:
		(cohort, cancer, platform, genome) = input.baseName.split('\\.')
		cutoff = CUTOFFS[platform]
		output = "${cohort}.${platform}.tsv.gz"
		if (cutoff)
			"""
			parse-variants --input ${input} --output ${output} \
				--genome ${genome.toLowerCase()} \
				--cutoff ${cutoff}
			"""
		else
			error "Invalid cutoff. Check platform: $platform"

}

// TODO rename lowercases
VARIANTS.into{ VARIANTS1; VARIANTS2; VARIANTS3; VARIANTS4; VARIANTS5 }


process formatSignature {
	tag "Prepare for signatures ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/signature", mode: "symlink", enabled: params.debug

	input:
		path input from VARIANTS1

	output:
		path output into VARIANTS_SIG

	script:
		(cohort, platform) = input.baseName.split('\\.')
		output = "${cohort}.${platform}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format signature
		"""

}

REGIONS_PREFIX = ['WXS': 'cds', 'WGS': 'wg']

process Signature {
	tag "Signatures ${cohort}"
	label "bgsignature"
	publishDir "${DEBUG_FOLDER}/signature", mode: "symlink", enabled: params.debug

	input:
		path input from VARIANTS_SIG

	output:
		tuple val(cohort), path(output) into SIGNATURES

	script:
		(cohort, platform) = input.baseName.split('\\.')
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

SIGNATURES.into{ SIGNATURES1; SIGNATURES2; SIGNATURES3; SIGNATURES4 }


process formatFML {
	tag "Prepare for FML ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/oncodrivefml", mode: "symlink", enabled: params.debug


	input:
		path input from VARIANTS2

	output:
		tuple val(cohort), path(output) into VARIANTS_FML

	script:
		(cohort, platform) = input.baseName.split('\\.')
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format fml
		"""

}

process OncodriveFML {
    tag "OncodriveFML ${cohort}"
    publishDir "${DEBUG_FOLDER}/oncodrivefml", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(input), path(signature)  from VARIANTS_FML.join(SIGNATURES1)
        path regions from REGIONS

    output:
        tuple val(cohort), path("out/*.tsv.gz") into OUT_ONCODRIVEFML

	script:
		// TODO add --debug in debug mode
		// TODO is the -c needed? We already have the
		"""
		oncodrivefml -i ${input} -e ${regions} --signature ${signature} \
			-c /oncodrivefml/oncodrivefml_v2.conf  --cores ${task.cpus} \
			-o out --debug
		"""
}


process formatCLUSTL {
	tag "Prepare for CLUSTL ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/oncodriveclustl", mode: "symlink", enabled: params.debug

	input:
		path input from VARIANTS3

	output:
		tuple val(cohort), path(output) into VARIANTS_CLUSTL

	script:
		(cohort, platform) = input.baseName.split('\\.')
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format clustl
		"""

}

process OncodriveCLUSTL {
    tag "OncodriveCLUSTL ${cohort}"
    publishDir "${DEBUG_FOLDER}/oncodriveclustl", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(input), path(signature)  from VARIANTS_CLUSTL.join(SIGNATURES2)
        path regions from REGIONS

    output:
        tuple val(cohort), path("out/elements_results.txt") into OUT_ONCODRIVECLUSTL
        tuple val(cohort), path("out/clusters_results.tsv") into OUT_ONCODRIVECLUSTL_CLUSTERS

	script:
		// TODO change kmer and other options for SKCM cancer type
		// TODO add --debug in debug mode
		"""
		oncodriveclustl -i ${input} -r ${regions} \
			-g hg38 -sim region_restricted -n 1000 -kmer 3 \
			-sig ${signature} --concatenate \
			-c ${task.cpus} \
			-o out
		"""
		// TODO is this needed ?
		//(cat {self.output_folder}/{self.name}/elements_results.txt | gzip > {self.out_file}) &&
        //    (cat {self.output_folder}/{self.name}/clusters_results.tsv | gzip > {self.output_folder}/{self.name}.clusters.gz) &&
        //    (cat {self.output_folder}/{self.name}/oncohortdrive_results.out | gzip > {self.output_folder}/{self.name}.oncohortdrive.gz)

}


process formatDNDSCV {
	tag "Prepare for DNDSCV ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/dndscv", mode: "symlink", enabled: params.debug

	input:
		path input from VARIANTS4

	output:
		tuple val(cohort), path(output) into VARIANTS_DNDSCV

	script:
		(cohort, platform) = input.baseName.split('\\.')
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format dndscv
		"""
}

process dNdScv {
    tag "dNdScv ${cohort}"
    publishDir "${DEBUG_FOLDER}/dndscv", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(input) from VARIANTS_DNDSCV

    output:
        tuple val(cohort), path("${cohort}.dndscv.tsv.gz") into OUT_DNDSCV

	script:
		"""
		Rscript /dndscv/dndscv.R \
			${input} ${cohort}.dndscv.tsv.gz \
			${cohort}.dndscv_annotmuts.tsv.gz \
			${cohort}.dndscv_genemuts.tsv.gz
		"""
}

OUT_DNDSCV.into{ OUT_DNDSCV1; OUT_DNDSCV2 }

process formatVEP {
	tag "Prepare for VEP ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/vep", mode: "symlink", enabled: params.debug

	input:
		path input from VARIANTS5

	output:
		tuple val(cohort), path(output) into VARIANTS_VEP

	script:
		(cohort, platform) = input.baseName.split('\\.')
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format vep
		"""

}

process VEP {
	tag "VEP ${cohort}"
	publishDir "${DEBUG_FOLDER}/vep", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(input) from VARIANTS_VEP

    output:
        tuple val(cohort), path(output) into OUT_VEP

	script:
		output = "${cohort}.vep.tsv.gz"
		"""
		vep -i ${input} -o STDOUT --assembly GRCh38 \
			--no_stats --cache --offline --symbol \
			--protein --tab --canonical \
			--dir ${params.datasets}/vep \
			| grep -v ^## | gzip > ${output}
		"""
}


process processVEPoutput {
	tag "Process vep output ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/vep", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(input) from OUT_VEP

    output:
        tuple val(cohort), path(output) into PARSED_VEP
        tuple val(cohort), path("${output}.stats.json") into STATS_VEP

	script:
		output = "${cohort}.tsv.gz"
		"""
		parse-vep --input ${input} --output ${output}
		"""
}


PARSED_VEP.into { PARSED_VEP1; PARSED_VEP2; PARSED_VEP3; PARSED_VEP4; PARSED_VEP5; PARSED_VEP6 }

process filterNonSynonymous {
	tag "Filter non synonymus ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/nonsynonymous", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(input) from PARSED_VEP1

    output:
        tuple val(cohort), path(output) into PARSED_VEP_NONSYNONYMOUS

	script:
		output = "${cohort}.vep_nonsyn.tsv.gz"
		"""
		parse-nonsynonymous --input ${input} --output ${output}
		"""
}


process formatSMRegions {
	tag "Prepare for SMRegions ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/smregions", mode: "symlink", enabled: params.debug

	input:
		tuple val(cohort), path(input) from PARSED_VEP_NONSYNONYMOUS

	output:
		tuple val(cohort), path(output) into VARIANTS_SMREGIONS

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format smregions
		"""
}

process SMRegions {
	tag "SMRegions ${cohort}"
	publishDir "${DEBUG_FOLDER}/smregions", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(input), path(signature)  from VARIANTS_SMREGIONS.join(SIGNATURES3)
        path regions from REGIONS

    output:
        tuple val(cohort), path(output) into OUT_SMREGIONS

	script:
		// TODO add --debug in debug mode
		output = "${cohort}.smregions.tsv.gz"
		"""
		smregions -m ${input} -e ${regions} \
			-r ${params.datasets}/smregions/regions_pfam.tsv \
			-s ${signature} --cores ${task.cpus} \
			-o ${output} --debug
		"""
}

OUT_SMREGIONS.into { OUT_SMREGIONS1; OUT_SMREGIONS2 }


process formatCBaSE {
	tag "Prepare for CBaSE ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/cbase", mode: "symlink", enabled: params.debug

	input:
		tuple val(cohort), path(input) from PARSED_VEP2

	output:
		tuple val(cohort), path(output) into VARIANTS_CBASE

	script:
		output = "${cohort}.in.tsv"
		"""
		format-variants --input ${input} --output ${output} \
			--format cbase
		"""
}

process CBaSE {
	tag "CBaSE ${cohort}"
	publishDir "${DEBUG_FOLDER}/cbase", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(input) from VARIANTS_CBASE

    output:
        tuple val(cohort), path(output) into OUT_CBASE

	script:
		output = "${cohort}.cbase.tsv.gz"
		"""
		python /cbase/cbase.py ${input} ${params.datasets}/cbase 0
		tail -n+2 q_values_output.txt | gzip > ${output}
		"""
}

process formatMutPanning {
	tag "Prepare for MutPanning ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/mutpanning", mode: "symlink", enabled: params.debug
	errorStrategy 'ignore' // TODO remove

	input:
		tuple val(cohort), path(input) from PARSED_VEP3

	output:
		tuple val(cohort), path(muts), path(samples) into VARIANTS_MUTPANNING

	script:
		muts = "${cohort}.in_muts.tsv"
		samples = "${cohort}.in_samples.tsv"
		"""
		format-variants --input ${input} --output ${muts} \
			--format mutpanning-mutations
		format-variants --input ${input} --output ${samples} \
			--format mutpanning-samples
		"""
}

process MutPanning {
	tag "MutPanning ${cohort}"
	publishDir "${DEBUG_FOLDER}/mutpanning", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(mutations), path(samples) from VARIANTS_MUTPANNING
        path regions from REGIONS

    output:
        tuple val(cohort), path("out/SignificanceFiltered/Significance${cohort}.txt") into OUT_MUTPANNING

	script:
		// TODO remove the creation of the out file or move to the container
		"""
		mkdir -p out/SignificanceFiltered
		echo "Name\tTargetSize\tTargetSizeSyn\tCount\tCountSyn\tSignificance\tFDR\n" \
			> out/SignificanceFiltered/Significance${cohort}.txt
		java -cp /mutpanning/MutPanning.jar MutPanning \
			out ${mutations} ${samples} ${params.datasets}/mutpanning/Hg19/
		"""
}


process formatHotMAPS {
	tag "Prepare for HotMAPS ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/hotmaps", mode: "symlink", enabled: params.debug

	input:
		tuple val(cohort), path(input) from PARSED_VEP4

	output:
		tuple val(cohort), path(output) into VARIANTS_HOTMAPS

	script:
		output = "${cohort}.in.maf"
		"""
		format-variants --input ${input} --output ${output} \
			--format hotmaps
		"""
}

process HotMAPS {
	tag "HotMAPS ${cohort}"
	publishDir "${DEBUG_FOLDER}/hotmaps", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(input), path(signatures) from VARIANTS_HOTMAPS.join(SIGNATURES4)
        path regions from REGIONS

    output:
        tuple val(cohort), path("*.out.gz") into OUT_HOTMAPS
        tuple val(cohort), path("*.clusters.gz") into OUT_HOTMAPS_CLUSTERS

	script:
		// TODO add --debug in debug mode
		"""
		/bin/sh /hotmaps/hotmaps.sh ${input} . ${signatures} \
			${params.datasets}/hotmaps ${task.cpus}
		"""
}


process Combination {
	tag "Combination ${cohort}"
	publishDir "${DEBUG_FOLDER}/combination", mode: "symlink", enabled: params.debug

    input:
        tuple val(cohort), path(fml), path(clustl), path(dndscv), path(smregions), path(cbase), path(mutpanning), path(hotmaps) from OUT_ONCODRIVEFML.join(OUT_ONCODRIVECLUSTL).join(OUT_DNDSCV1).join(OUT_SMREGIONS1).join(OUT_CBASE).join(OUT_MUTPANNING).join(OUT_HOTMAPS)

    output:
        tuple val(cohort), path("${cohort}.05.out.gz") into OUT_COMBINATION

	script:
		// TODO add --debug in debug mode
		"""
		intogen-combine -o ${cohort} \
			--oncodrivefml ${fml} \
			--oncodriveclustl ${clustl} \
			--dndscv ${dndscv} \
			--smregions ${smregions} \
			--cbase ${cbase} \
			--mutpanning ${mutpanning} \
			--hotmaps ${hotmaps}
		"""

}


process formatdeconstructSigs {
	tag "Prepare for deconstructSigs ${cohort}"
	label "core"
	publishDir "${DEBUG_FOLDER}/deconstructSigs", mode: "symlink", enabled: params.debug

	input:
		tuple val(cohort), path(input) from PARSED_VEP5

	output:
		tuple val(cohort), path(output) into VARIANTS_DECONSTRUCTSIGS

	script:
		output = "${cohort}.in.tsv.gz"
		"""
		format-variants --input ${input} --output ${output} \
			--format deconstructsigs
		"""
}

VARIANTS_DECONSTRUCTSIGS.into{ VARIANTS_DECONSTRUCTSIGS1; VARIANTS_DECONSTRUCTSIGS2 }

process deconstructSigs {
	tag "deconstructSigs ${cohort}"
	publishDir "${DEBUG_FOLDER}/deconstructSigs", mode: "symlink", enabled: params.debug

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


// TODO add a process to combine all stats

process CohortCounts {
	tag "Count variants ${cohort}"

    input:
        path(input) from PARSED_INPUT2

    output:
		tuple val(cohort), path("*.counts") into COHORT_COUNTS

	script:
		// TODO add cancer type
		(cohort, platform, genome) = input.baseName.split('\\.')
		"""
		variants=`zcat ${input} |tail -n+2 |wc -l`
		#C=1; for i in \$(head ${input} -n 1) ; do if [ \$i == "SAMPLE" ] ; then break ; else C=\$(( \$C + 1 )) ; fi ; done
		samples=`zcat ${input} |tail -n+2 |cut -f1 |uniq |sort -u| wc -l`
		echo "${cohort}\t${platform}\t\$variants\t\$samples" > ${cohort}.counts
		"""
}

COHORT_COUNTS.into{ COHORT_COUNTS1; COHORT_COUNTS2 }
COHORT_COUNTS_LIST = COHORT_COUNTS1.map{ it -> it[1] }

process CohortSummary {
	tag "Count variants ${cohort}"
	publishDir "${OUTPUT}", mode: "copy"

    input:
        path(input) from COHORT_COUNTS_LIST.collect()

    output:
		path(output) into COHORT_SUMMARY

	script:
		output="cohorts.tsv"
		"""
		echo 'COHORT\tPLATFORM\tMUTATIONS\tSAMPLES' > ${output}
		cat ${input} >> ${output}
		"""
}


MUTATIONS_INPUTS = PARSED_VEP6.map { it -> it[1] }

process MutationsSummary {
	tag "Mutations"
	publishDir "${OUTPUT}", mode: "copy"
	label "core"

    input:
        path(input) from MUTATIONS_INPUTS.collect()

    output:
		path(output) into MUTATIONS_SUMMARY

	script:
		output="mutations.tsv"
		"""
		mutations-summary --output ${output} \
			${input}
		"""
}
//FIXME remove 2 from name
process DriverDiscovery2 {
	tag "Driver discovery ${cohort}"
	publishDir "${DEBUG_FOLDER}/drivers", mode: "symlink", enabled: params.debug
	label "core"

    input:
        tuple val(cohort), path(combination), path(deconstruct_in), path(sig_likelihood), path(smregions), path(clustl_clusters), path(hotmaps_clusters), path(dndscv) from OUT_COMBINATION.join(VARIANTS_DECONSTRUCTSIGS2).join(OUT_DECONSTRUCTSIGS_SIGLIKELIHOOD).join(OUT_SMREGIONS2).join(OUT_ONCODRIVECLUSTL_CLUSTERS).join(OUT_HOTMAPS_CLUSTERS).join(OUT_DNDSCV2)

    output:
		path(output) into DRIVERS

	script:
		output = "${cohort}.drivers.tsv"
		// FIXME get tumor type
		"""
		drivers-discovery --output ${output} \
			--combination ${combination} \
			--mutations ${deconstruct_in} \
			--sig_likelihood ${sig_likelihood} \
			--smregions ${smregions} \
			--clustl_clusters ${clustl_clusters} \
			--hotmaps ${hotmaps_clusters} \
			--dndscv ${dndscv} \
			--ctype ACC --cohort ${cohort}
		"""
}
//FIXME remove 2 from name
process DriverSummary2 {
	tag "Driver summary"
	publishDir "${OUTPUT}", mode: "copy"
	label "core"

    input:
        path (input) from DRIVERS.collect()
        path (mutations) from MUTATIONS_SUMMARY
        path (cohort) from COHORT_SUMMARY

    output:
		path("*.tsv") into DRIVERS_SUMMARY

	script:
		// FIXME get tumor type
		"""
		drivers-summary \
			--mutations ${mutations} \
			--cohorts ${cohort} \
			${input}
		"""
}
