//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FORMAT            as FORMAT_MUTPANNING_SAMPLE    } from '../../../modules/local/core/format_variants/main'
include { FORMAT            as FORMAT_MUTPANNING_MUTS      } from '../../../modules/local/core/format_variants/main'
include { MUTPANNING                                       } from '../../../modules/local/mutpanning/main'


workflow MUTPANNING {
    take:
    parsed_vep
    
    main:
    ch_versions     = Channel.empty()
    format_sample   = 'mutpanning-sample'
    format_muts     = 'mutpanning-mutations'

    FORMAT_MUTPANNING_SAMPLE(
        parsed_vep, 
        format_sample
    )

    FORMAT_MUTPANNING_MUTS(
        parsed_vep,
        format_muts
    )

    FORMAT_MUTPANNING_MUTS.out.format
                                .set{ mutations }
    FORMAT_MUTPANNING_SAMPLE.out.format
                                .set{ samples }

    MUTPANNING(
        mutations,
        samples
    )

    ch_versions = ch_versions.mix(FORMAT_MUTPANNING_MUTS.out.versions)
    ch_versions = ch_versions.mix(FORMAT_MUTPANNING_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(MUTPANNING.out.versions)

    emit:
    results     = MUTPANNING.out.out_mutpanning
    versions    = ch_versions                // channel: [ versions.yml ]
}

