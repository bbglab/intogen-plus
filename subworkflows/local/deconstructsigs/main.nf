//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FORMAT               as FORMAT_DECONSTRUCSIGS       } from '../../../modules/local/core/format_variants/main'
include { DECONSTRUCSIGS                                      } from '../../../modules/local/deconstructsigs/main'


workflow DECONSTRUCSIGS {
    take:
    parsed_vep
    
    main:
    ch_versions = Channel.empty()
    format      = 'deconstructsigs'

    FORMAT_DECONSTRUCSIGS(
        parsed_vep,
        format
    )

    DECONSTRUCSIGS(
        FORMAT_DECONSTRUCSIGS.out.format
    )

    ch_versions = ch_versions.mix(FORMAT_DECONSTRUCSIGS.out.versions)
    ch_versions = ch_versions.mix(DECONSTRUCSIGS.out.versions)

    emit:
    results             = DECONSTRUCSIGS.out.out_DECONSTRUCSIGS
    likelihood_results  = DECONSTRUCSIGS.out.out_DECONSTRUCSIGS_likelihood
    versions            = ch_versions                // channel: [ versions.yml ]
}