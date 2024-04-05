//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FORMAT        as FORMAT_HOTMAPS       } from '../../../modules/local/core/format_variants/main'
include { HOTMAPS       as HOTMAPS              } from '../../../modules/local/hotmaps/main'


workflow HOTMAPS {
    take:
    parsed_vep

    main:
    ch_versions = Channel.empty()
    format      = 'hotmaps'

    FORMAT_HOTMAPS(
        parse_vep,
        format
    )

    HOTMAPS(
        profile,
        FORMAT_HOTMAPS.out.format
    )

    ch_versions = ch_versions.mix(FORMAT_HOTMAPS.out.versions)
    ch_versions = ch_versions.mix(HOTMAPS.out.versions)

    emit:
    results             = HOTMAPS.out.out_hotmaps
    clusters_results    = HOTMAPS.out.out_hotmaps_clusters
    versions            = ch_versions                // channel: [ versions.yml ]
}