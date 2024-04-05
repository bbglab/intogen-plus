//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { GET_FIELD         as GET_TUMORTYPE        } from '../../../modules/local/get_field/main'
include { GET_FIELD         as GET_PLATFORM         } from '../../../modules/local/get_field/main'
include { GET_FIELD         as GET_GENOME           } from '../../../modules/local/get_field/main'

include { PROCESS_VARIANTS                          } from '../../../modules/local/core/process_variants/main'
include { COUNTS_VARIANTS                           } from '../../../modules/local/countVariants/main'


workflow PARSE_INPUT {
    take:
    cohort
    
    main:
    ch_versions      = Channel.empty()
    
    fieldGenome      = 'GENOMEREF'
    fieldPlatform    = 'PLATFORM'
    fieldTumortype   = 'CANCER'


    GET_TUMORTYPE( cohort, fieldTumortype )
    GET_GENOME( cohort, fieldGenome )
    GET_PLATFORM ( cohort, fieldPlatform )


    PROCESS_VARIANTS( 
        cohort,
        GET_PLATFORM.out.platform,
        GET_GENOME.out.genome
    )

    COUNTS_VARIANTS(
        cohort,
        GET_TUMORTYPE.out.tumor_type,
        GET_PLATFORM.out.platform,
    )


    ch_versions = ch_versions.mix(FORMAT_MUTPANNING_MUTS.out.versions)
    ch_versions = ch_versions.mix(FORMAT_MUTPANNING_SAMPLE.out.versions)
    ch_versions = ch_versions.mix(MUTPANNING.out.versions)

    emit:
    variants        = PROCESS_VARIANTS.out.out_variants
    tumor_type      = GET_TUMORTYPE.out.out_field
    platform        = GET_PLATFORM.out.out_field
    genome          = GET_GENOME.out.out_field
    variants_count  = COUNTS_VARIANTS.out.out_count
    versions        = ch_versions                // channel: [ versions.yml ]
}

