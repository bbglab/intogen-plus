//
// Subworkflow with functionality specific to the bbglab/intogen-plus pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FORMAT        as FORMAT_DNDSCV        } from '../../../modules/local/core/format_variants/main'
include { DNDSCV                                } from '../../../modules/local/dndscv/main'


workflow DNDSCV {
    take:
    variants
    
    main:
    ch_versions = Channel.empty()
    format      = 'dndscv'

    FORMAT_DNDSCV(
        variants, 
        format
    )

    DNDSCV(
        FORMAT_DNDSCV.out.format
    )

    ch_versions = ch_versions.mix(FORMAT_DNDSCV.out.versions)
    ch_versions = ch_versions.mix(DNDSCV.out.versions)

    emit:
    results             = DNDSCV.out.out_dndscv
    annot_muts_results  = DNDSCV.out.out_dndscv_annotated_mutations
    genemuts_results    = DNDSCV.out.out_dndscv_gene_mutations
    versions            = ch_versions                // channel: [ versions.yml ]
}

