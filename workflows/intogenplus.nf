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

include { PARSE_INPUT            as PARSE             } from '../subworkflow/local/parseInput/main'
include { PREPROCESS_VARIANTS    as PREPROCESS        } from '../subworkflow/local/preprocess/main'

include { GET_PROFILE            as PROFILE           } from '../subworkflow/local/get_profile/main'
include { VEP                    as VEP               } from '../subworkflow/local/vep/main' // check if there's vep nf-core

include { ONCODRIVEFML                                } from '../subworkflow/local/oncodrivefml/main'
include { ONCODRIVECLUSTL                             } from '../subworkflow/local/oncodriveclustl/main'
include { SMREGIONS                                   } from '../subworkflow/local/smregions/main'
include { DNDSCV                                      } from '../subworkflow/local/dndscv/main'
include { CBASE                                       } from '../subworkflow/local/cbase/main'
include { MUTPANNING                                  } from '../subworkflow/local/mutpanning/main'
include { HOTMAPS                                     } from '../subworkflow/local/hotmaps/main'
// include { ONCODRIVE3D                                 } from '../subworkflow/local/oncodrive3d/main'
// include { SEISMIC                                     } from '../subworkflow/local/seismic/main'

include { DECONSTRUCSIGS                              } from '../subworkflow/local/deconstructsigs/main'

include { INTOGENCOMBINATION     as COMBINATION       } from '../module/local/combination/main'
include { CORE_DRIVER_DISCOVERY  as DISCOVERY         } from '../module/local/core/driver_discovery/main'

include { INTOGENSUMMARY         as ALLSTATS           } from '../subworkflow/local/intogen_summary/main'

include { BOOSTDM_CONNECTION     as BOOSTDM_CONN      } from '../subworkflow/local/intogen_boostdm_connection/main' 


// include { paramsSummaryMap                            } from 'plugin/nf-validation'
// include { paramsSummaryMultiqc                        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                      } from '../subworkflows/local/utils_nfcore_intogenplus_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INTOGENPLUS {

    take:

    input
    regions

    main:

    ch_versions = Channel.empty()

    //
    // STEP0: SUBWORKFLOW: Parsing the input
    //
    PARSE(input)    
    PARSE.out.cohorts
                .flatten()
                .map{ it -> [it.baseName.split('\\.')[0], it]}
                .set{ cohort }

    //
    // STEP1: Preprocess
    //
    PREPROCESS(cohort)
    PREPROCESS.out.variants
                    .set { varsXcohort }
    PREPROCESS.out.tumor_type
                    .set { ttype }
    PREPROCESS.out.platform
                    .set { platform }
    PREPROCESS.out.variants_count
                    .map{ it -> it[1] }
                    .set { variantsCount }

    //
    // STEP1.1: Preprocess: Get Profile
    //
    PROFILE( varsXcohort )
    PROFILE.out.profile.set{ profile }

    //
    // STEP1.2: Preprocess: Annotate Variants
    //
    VEP( varsXcohort )
    VEP.out.parsed_vep
                .set { varsVepped }
    
    DECONSTRUCSIGS( varsVepped )

    //
    // STEP2: Run methods
    //
    DNDSCV( varsXcohort )
    ONCODRIVEFML( varsXcohort, profile, regions )
    ONCODRIVECLUSTL( varsXcohort, profile, regions, ttype)
    CBASE( varsVepped )
    MUTPANNING( varsVepped )
    HOTMAPS( varsVepped )    
    SMREGIONS( varsVepped, profile, regions )

    //
    // STEP3: MODULE: Run combination
    //
    COMBINATION(
        ONCODRIVEFML.out.results,
        ONCODRIVECLUSTL.out.results,
        DNDSCV.out.results,
        SMREGIONS.out.results,
        CBASE.out.results,
        MUTPANNING.out.results,
        HOTMAPS.out.results,
    )

    //
    // STEP4: Discover Drivers and Summary
    //
    DISCOVERY(
        COMBINATION.out.results,
        DECONSTRUCSIGS.out.results,
        DECONSTRUCSIGS.out.likelihood_results,
        SMREGIONS.out.results,
        ONCODRIVECLUSTL.out.clusters_results,
        HOTMAPS.out.clusters_results,
        DNDSCV.out.results,
        tumor_type
    )

    DISCOVERY.out.drivers
                    .collect()
                    .set{ allDrivers }
    DISCOVERY.out.drivers   
                    .collect()
                    .set{ allDriversVet }
    varsVepped
        .collect
        .set{ allVariants }

    variantsCount
        .collect
        .set{ allVarsCount}
    
    
    ALLSTATS(
        allDrivers,
        allDriversVet,
        allVariants,
        allVarsCount,
    )

    //
    // STEP5: Generate BoostDM connection deps
    //
    BOOSTDM_CONN(
        allVariants,
        ALLSTATS.out.summaryDrivers
    )
    
    

    // collect ch_versions
    ch_versions = ch_versions.mix(PARSE.out.versions)
    ch_versions = ch_versions.mix(PREPROCESS.out.versions)
    ch_versions = ch_versions.mix(GET_PROFILE.out.versions)
    ch_versions = ch_versions.mix(VEP.out.versions)
    ch_versions = ch_versions.mix(ONCODRIVECLUSTL.out.versions)
    ch_versions = ch_versions.mix(ONCODRIVEFML.out.versions)
    ch_versions = ch_versions.mix(SMREGIONS.out.versions)
    ch_versions = ch_versions.mix(DNDSCV.out.versions)
    ch_versions = ch_versions.mix(CBASE.out.versions)
    ch_versions = ch_versions.mix(MUTPANNING.out.versions)
    ch_versions = ch_versions.mix(HOTMAPS.out.versions)
    ch_versions = ch_versions.mix(DECONSTRUCSIGS.out.versions)
    ch_versions = ch_versions.mix(COMBINATION.out.versions)
    ch_versions = ch_versions.mix(DISCOVERY.out.versions)
    ch_versions = ch_versions.mix(ALLSTATS.out.versions)
    ch_versions = ch_versions.mix(BOOSTDM_CONN.out.versions)


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowDeepcsa.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowDeepcsa.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()


    emit:
    multiqc_report = 4.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
