//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FILTER_NONSYNONYMOUS  as FILTER_NONSYN        } from '../../../modules/local/core/parse_nonsynonymous/main'
include { FORMAT                as FORMAT_SMREGIONS     } from '../../../modules/local/core/format_variants/main'
include { SMREGIONS                                     } from '../../../modules/local/smregions/main'


workflow SMREGIONS {
    take:
    parsed_vep
    profile
    regions
    
    main:
    ch_versions = Channel.empty()
    format      = 'smregions'

    FILTER_NONSYN(
        parsed_vep,
    )

    FORMAT_SMREGIONS(
        FILTER_NONSYN.out.parsed_vep_nonsyn, 
        format
    )

    SMREGIONS(
        profile
        FORMAT_SMREGIONS.out.format
        )

    ch_versions = ch_versions.mix(FILTER_NONSYN.out.versions)
    ch_versions = ch_versions.mix(FORMAT_SMREGIONS.out.versions)
    ch_versions = ch_versions.mix(SMREGIONS.out.versions)

    emit:
    results     = SMREGIONS.out.out_SMREGIONS
    versions    = ch_versions                // channel: [ versions.yml ]
}