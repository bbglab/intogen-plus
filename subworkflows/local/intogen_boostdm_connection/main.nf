//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FILTERMNVS                                 } from '../../../modules/local/core/parse_mnvs/main'
include { DRIVERSSATURATION        as SATURATION     } from '../../../modules/local/core/drivers_saturation/main'


workflow PARSE_INPUT {
    take:
    all_variants
    drivers_summary

    main:
    ch_versions = Channel.empty()

    FILTERMNVS( all_variants )

    SATURATION( driver_summary )

    ch_versions = ch_versions.mix(FILTERMNVS.out.versions)
    ch_versions = ch_versions.mix(SATURATION.out.versions)
    
    emit:
    mnvsFilter = FILTERMNVS.out.mnvs_filterç
    driversSaturation = SATURATION.out.drivers_saturation
    versions = ch_versions                // channel: [ versions.yml ]
}

