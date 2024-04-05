//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FORMAT        as FORMAT_CBASE     } from '../../../modules/local/core/format_variants/main'
include { CBASE                             } from '../../../modules/local/cbase/main'


workflow CBASE {
    take:
    parsed_vep
    
    main:
    ch_versions = Channel.empty()
    format      = 'cbase'

    FORMAT_CBASE(
        parsed_vep, 
        format
    )

    CBASE(
        FORMAT_CBASE.out.format
    )

    ch_versions = ch_versions.mix(FORMAT_CBASE.out.versions)
    ch_versions = ch_versions.mix(CBASE.out.versions)

    emit:
    results     = CBASE.out.out_cbase
    versions    = ch_versions                // channel: [ versions.yml ]
}

