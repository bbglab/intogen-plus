/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowDeepcsa.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PARSE_INPUT            as PARSE_INPUT       } from '../subworkflow/local/parseInput/main'
include { PREPROCESS_VARIANTS    as PREPROCESS        } from '../subworkflow/local/preprocess_variants/main'

include { GET_PROFILE            as GET_PROFILE       } from '../subworkflow/local/get_profile/main'
include { VEP                    as VEP               } from '../subworkflow/local/vep/main' // check if there's vep nf-core

include { ONCODRIVEFML           as ONCODRIVEFML      } from '../subworkflow/local/oncodrivefml/main'
include { ONCODRIVECLUSTL        as ONCODRIVECLUSTL   } from '../subworkflow/local/oncodriveclustl/main'
include { SMREGIONS              as SMREGIONS         } from '../subworkflow/local/smregions/main'
include { DNDSCV                 as DNDSCV            } from '../subworkflow/local/dndscv/main'
include { CBASE                  as CBASE             } from '../subworkflow/local/cbase/main'
include { MUTPANNING             as MUTPANNING        } from '../subworkflow/local/mutpanning/main'
include { HOTMAPS                as HOTMAPS           } from '../subworkflow/local/hotmaps/main'
// include { ONCODRIVE3D           as ONCODRIVE3D        } from '../subworkflow/local/oncodrive3d/main'
// include { SEISMIC               as SEISMIC            } from '../subworkflow/local/seismic/main'

include { DECONSTRUCSIGS         as DECONSTRUCSIGS    } from '../subworkflow/local/deconstructsigs/main'

include { INTOGENCOMBINATION     as COMBINATION       } from '../module/local/combination/main'
include { CORE_DRIVER_DISCOVERY  as DISCOVERY         } from '../module/local/core/driver_discovery/main'

include { INTOGENSUMMARY         as SUMMARY           } from '../subworkflow/local/intogen_summary/main'

include { BOOSTDM_CONNECTION     as BOOSTDM_CONN      } from '../subworkflow/local/intogen_boostdm_connection/main' 


include { paramsSummaryMap                            } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                      } from '../subworkflows/local/utils_nfcore_intogenplus_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INTOGENPLUS {

    take:

    INPUT    = Channel.fromPath(params.input.tokenize())
    REGIONS  = Channel.value("${params.datasets}/regions/cds.regions.gz")

    main:

    ch_versions = Channel.empty()

    //
    // STEP0: SUBWORKFLOW: Parsing the input
    //
    PARSE(INPUT)    

    PARSE.out.cohorts
                .flatten()
                .map{ it -> [it.baseName.split('\\.')[0], it]}
                .set{ cohort }

    //
    // STEP1: Preprocess and annotate input
    //
    PREPROCESS_VARIANTS(cohort)
    
    PREPROCESS_VARIANTS.out.variants
                                .set { varsXcohort }
    PREPROCESS_VARIANTS.out.tumor_type
                                .set { ttype }
    PREPROCESS_VARIANTS.out.platform
                                .set { platform }


    GET_PROFILE( varsXcohort )
    GET_PROFILE.out.profile.set{ profile }


    VEP( varsXcohort )
    VEP.out.parsed_vep
                .set { varsVepped }
    
    //
    // STEP2: Run methods
    //
    DNDSCV( varsXcohort )

    ONCODRIVEFML( varsXcohort, profile, REGIONS )
    
    ONCODRIVECLUSTL( varsXcohort, profile, REGIONS, ttype)
    
    CBASE( varsVepped )

    MUTPANNING( varsVepped )

    HOTMAPS( varsVepped )    
    
    DECONSTRUCSIGS( varsVepped )

    SMREGIONS( varsVepped, profile, REGIONS )
    
    //
    // STEP3: MODULE: Run combination
    //
    COMBINATION()

    //
    // STEP4: Discover Drivers and do summary
    //
    DISCOVERY()

    SUMMARY()

    //
    // STEP5: Produce file for BoostDM connection
    //
    BOOSTDM_CONN()
    
    

    // collect ch_versions
    ch_versions = ch_versions.mix(PARSE.out.versions)
    ch_versions = ch_versions.mix(PREPROCESS_VARIANTS.out.versions)
    ch_versions = ch_versions.mix(GET_PROFILE.out.versions)
    ch_versions = ch_versions.mix(VEP.out.versions)
    ch_versions = ch_versions.mix(ONCODRIVECLUSTL.out.versions)
    ch_versions = ch_versions.mix(ONCODRIVEFML.out.versions)
    ch_versions = ch_versions.mix(SMREGIONS.out.versions)
    ch_versions = ch_versions.mix(DNDSCV.out.versions)
    ch_versions = ch_versions.mix(CBASE.out.versions)
    ch_versions = ch_versions.mix(MUTPANNING.out.versions)
    ch_versions = ch_versions.mix(HOTMAPS.out.versions)
    ch_versions = ch_versions.mix(HOTMAPS.out.versions)
    ch_versions = ch_versions.mix(COMBINATION.out.versions)
    ch_versions = ch_versions.mix(DISCOVERY.out.versions)
    ch_versions = ch_versions.mix(SUMMARY.out.versions)
    ch_versions = ch_versions.mix(BOOSTDM_CONN.out.versions)
  
    


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
