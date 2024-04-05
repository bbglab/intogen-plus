//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FORMAT                as FORMAT_CLUSTL        } from '../../../modules/local/core/format_variants/main'
include { ONCODRIVECLUSTL       as CLUSTL               } from '../../../modules/local/oncodriveclustl/main'


workflow ONCODRIVECLUSTL {
    take:
    tumor_type
    variants
    profile
    regions

    
    main:
    ch_versions = Channel.empty()
    format      = 'clustl'

    FORMAT_CLUSTL(
        variants, 
        format
    )

    CLUSTL(
        regions,
        tumor_type,
        profile
        FORMAT_CLUST.out.format
    )

    ch_versions = ch_versions.mix(FORMAT_CLUST.out.versions)
    ch_versions = ch_versions.mix(CLUSTL.out.versions)

    emit:
    results             = CLUSTL.out.out_clustl
    clusters_results    = CLUSTL.out.out_clusters_clustl
    versions            = ch_versions                // channel: [ versions.yml ]
}

