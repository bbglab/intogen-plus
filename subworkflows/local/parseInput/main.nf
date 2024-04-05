//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { OPENVARIANT_CAT           as CAT                  } from '../../../modules/local/openvariant/cat/main'
include { OPENVARIANT_GROUPBY       as GROUPBY              } from '../../../modules/local/openvariant/groupby/main'


workflow PARSE_INPUT {
    take:
    input

    main:
    ch_versions = Channel.empty()

    if ( input.toRealPath().toFile().isDirectory() || input.endsWith(".bginfo" )) {
        
        GROUPBY(input)
        GROUPBY.out.parsed_output
                        .set{ cohort }

        ch_versions = ch_versions.mix(GROUPBY.out.versions)
    }
    else {
        
        CAT(input)
        CAT.out.parsed_output
                        .set{ cohort }
        
        ch_versions = ch_versions.mix(CAT.out.versions)
    }

    emit:
    cohorts     = cohorts
    versions    = ch_versions                // channel: [ versions.yml ]
}

