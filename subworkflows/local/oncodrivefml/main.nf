//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FORMAT            as FORMAT_FML       } from '../../../modules/local/core/format_variants/main'
include { ONCODRIVEFML      as FML              } from '../../../modules/local/oncodrivefml/main'


workflow ONCODRIVEFML {
    take:
    variants
    profile
    regions
    
    main:
    ch_versions = Channel.empty()
    format      = 'fml'

    FORMAT_FML(
        variants, 
        format
    )

    FML(
        profile
        FORMAT_FML.out.format,
        regions
    )

    ch_versions = ch_versions.mix(FORMAT_FML.out.versions)
    ch_versions = ch_versions.mix(FML.out.versions)

    emit:
    results     = FML.out.out_FML
    versions    = ch_versions                // channel: [ versions.yml ]
}

