//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FORMAT                as FORMAT_PROFILE       } from '../../../modules/local/core/format_variants/main'
include { COMPUTE_PROFILE       as COMPUTE_PROFILE      } from '../../../modules/local/bgsignature/main'
include { PARSE_PROFILE         as PARSE_PROFILE        } from '../../../modules/local/core/parse_profile/main'


workflow GET_PROFILE {
    take:
    variants

    main:
    ch_versions = Channel.empty()


    FORMAT_PROFILE( 
        variants 
    )
    
    COMPUTE_PROFILE( 
        FORMAT_PROFILE.out.variantsProfile 
    )
    
    PARSE_PROFILE( 
        COMPUTE_PROFILE.out.profile 
    )

    ch_versions = ch_versions.mix(FORMAT_PROFILE.out.versions)
    ch_versions = ch_versions.mix(COMPUTE_PROFILE.out.versions)
    ch_versions = ch_versions.mix(PARSE_PROFILE.out.versions)

    emit:
    profile             = COMPUTE_PROFILE.out.profile
    boostdmProfile      = PARSE_PROFILE.out.boostdm_profile
    versions            = ch_versions                // channel: [ versions.yml ]
}

