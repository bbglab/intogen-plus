//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FORMAT                    as FORMAT_VEP       } from '../../../modules/local/core/format_variants/main'
include { VEP                       as VEP              } from '../../../modules/local/vep/main'
include { PROCESS_VEP_OUTPUT        as POSTPROCESS      } from '../../../modules/local/core/parse_vep/main'


workflow VEP {
    take:
    variants
    
    main:
    ch_versions = Channel.empty()
    format      = 'vep'

    FORMAT_VEP(
        variants,
        format
    )

    VEP(
        FORMAT_VEP.out.format, 
    )

    POSTPROCESS(
        VEP.out.out_vep
        )

    ch_versions = ch_versions.mix(FORMAT_VEP.out.versions)
    ch_versions = ch_versions.mix(VEP.out.versions)
    ch_versions = ch_versions.mix(POSTPROCESS.out.versions)

    emit:
    parsed_vep     = POSTPROCESS.out.out_vep_postprocess
    versions    = ch_versions                // channel: [ versions.yml ]
}