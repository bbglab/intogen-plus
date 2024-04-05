//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MUTATION_SUMMARY } from '../../../modules/local/core/mutation_summary/main'
include { COHORT_SUMMARY   } from '../../../modules/local/core/cohort_summary/main'
include { DRIVERS_SUMMARY  } from '../../../modules/local/core/drivers_summary/main'


workflow INTOGENSUMMARY {
    take:
    all_drivers
    all_drivers_vet
    all_variants
    all_variants_count

    main:
    ch_versions = Channel.empty()


    COHORT_SUMMARY( 
        all_variants_count 
    )
    
    MUTATION_SUMMARY( 
        all_variants 
    )
    
    DRIVERS_SUMMARY(
        all_drivers,
        all_drivers_vet,
        MUTATION_SUMMARY.out.mutsummary,
        COHORT_SUMMARY.out.cohortsummary
    )

    ch_versions = ch_versions.mix(COHORT_SUMMARY.out.versions)
    ch_versions = ch_versions.mix(MUTATION_SUMMARY.out.versions)
    ch_versions = ch_versions.mix(DRIVERS_SUMMARY.out.versions)

    emit:
    summaryDrivers      = DRIVERS_SUMMARY.out.drivers_summary
    uniqueDrivers       = DRIVERS_SUMMARY.out.drivers_unique
    unfilteredDrivers   = DRIVERS_SUMMARY.out.drivers_unfiltered 
    versions = ch_versions                // channel: [ versions.yml ]
}

