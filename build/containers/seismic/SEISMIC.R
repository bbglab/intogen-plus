#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(fitdistrplus))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plyranges))
suppressPackageStartupMessages(library(yaml))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(foreach))
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(stringi))
suppressPackageStartupMessages(library(inline))
suppressPackageStartupMessages(library(Rcpp))
suppressPackageStartupMessages(library(cowplot))

main <- function(){
  args <- commandArgs(trailingOnly = TRUE)

  config <- read_yaml("/seismic/config/config.yaml")
  
  config$mutations_path <- args[1]
  config$out_base_name <- str_split_i(config$mutations_path, "\\.", 1)
  config$cancer_type <- args[2]
  config$seq_type <- unlist(stri_split_fixed(config$out_base_name, "_", 3))[2]
  
  config$reference_genome <- args[3]
  no_of_cores <- as.integer(args[4])
  
  #print_logo()
  
  # If we are not dealing with CDSs, all mutation effects will be 'na', so change effects_to_keep to that to not filter anything out.
  if(!config$annotate_cds_effects) config$effects_to_keep <- 'na'
  
  genome <- load_bsgenome(config$reference_genome)
  fig_dir <- paste0(config$out_dir, 'figures/')
  system(paste("mkdir -p", config$out_dir, fig_dir))
  
  cat("Loading regions\n")
  # Read all regions that should be tested, e.g. cds regions. Could also be promoters, etc.
  test_regions.gr <- read_test_regions(config$test_regions_path, config$reference_genome)
  
  # Start a list of test regions (typically genes) with everything that might be analysed based on every region that has annotations
  # Filter to remove undesired regions later
  test_region_list <- test_regions.gr$test_region %>% unique()
  
  # Load region annotations, or annotate if there is no pre-existing file.
  mut_effects.df <- get_mut_effects(test_regions.gr, genome, config$reference_genome, config$test_regions_path, config$annotate_cds_effects, no_of_cores)
  
  cat("Loading mutations\n")
  # Read mutations. Use cancer_column_name and patient_column_name arguments if specified in config file. If not, leave blank to use function defaults.
  all_mutations.gr <- do.call(read_mutations,
                              list(path = config$mutations_path,
                                   cancer_type = config$cancer_type,
                                   genome = genome) %>%
                                {if(!is.null('config$maf_cancer_column_name')){c(., cancer_colname = config$maf_cancer_column_name)} else .} %>%
                                {if(!is.null('config$maf_patient_column_name')){c(., patient_colname = config$maf_patient_column_name)} else .}
  )
  
  no_of_patients <- all_mutations.gr$sampleID %>% unique() %>% length()
  
  cat("Annotating mutations\n")
  annotated_mutations.gr <- annotate_mutations(all_mutations.gr, mut_effects.df)
  effect_filtered_mutations.gr <- filter_mutation_effect_keep(annotated_mutations.gr, config$effects_to_keep)
  mutations_with_effects_for_plotting.df <- annotated_mutations.gr %>% as_tibble() %>% select(sampleID, base_no, codon_no, test_region, effect, refnuc, varnuc, strand)
  
  
  
  
  cat("Calculating mutational frequencies\n")
  if(config$sig_type == "patient"){
    #if(config$seq_type == "WXS"){
    #  cat("Patient-specific signatures not available for WXS sequencing data. Change sig_type accordingly\n")
    #  quit()
    #}
    sig_string <- "patient_sig_"
  } else if(config$sig_type == "cohort"){
    sig_string <- "cohort_sig_"
  } else if(config$sig_type == "flat"){
    sig_string <- "flat_sig_"
  } else {
    cat("wrong/no sig_type specified\n")
    quit()
  }
  mutfreqs.df <- calculate_mutfreqs(all_mutations.gr, config$sig_type, config$seq_type, genome, config$reference_genome, config$wxs_covered_regions_path)
  combined_mutfreqs.df <- combine_mutfreq_varnucs(mutfreqs.df)
  
  
  
  ################################################
  # Scaling expectations preparation
  ################################################
  #
  # Scaling! Currently, the scaling types are seq_type-specific
  scale_with_regression <- FALSE
  scale_silent_muts <- FALSE
  # Scaling options for WXS - linear regressions with replication timing, expression data, or both
  if(config$seq_type == "WXS"){
    if(config$scale_exp_to_repl_reg || config$scale_exp_to_expr_reg){
      cat("Performing regression for scaling expected mutations\n")
      scale_with_regression <- TRUE
      reglist <- initiate_regression_list(test_regions.gr, all_mutations.gr)
      
      if(config$scale_exp_to_repl_reg){
        cat("  Replication timing\n")
        reglist <- reglist %>%  scale_exp_to_repl_reg(config$repl_timing_path, test_regions.gr)
      }
      if(config$scale_exp_to_expr_reg){
        cat("  Expression\n")
        reglist <- reglist %>% scale_exp_to_expr_reg(config$expression_path, config$cancer_type)
        if(!("log_expr" %in% reglist$variables)){
          cat("    Expression scaling skipped\n")
        }
      }
      
      reglist <- make_scaling_regression(reglist, config$min_test_region_length_for_regression)
      scaling_string <- reglist$scale_with_regression_string
      
      # Filter out test regions that we don't have regression data for
      test_region_list <- test_region_list[test_region_list %in% reglist$test_region_mutrate$test_region]
    }
  } else if (config$seq_type == "WGS"){ # Scaling options for WGS - only silent mutations, currently.
    if(config$scale_exp_with_silent_muts){
      cat("Analysing silent mutations for scaling expected mutations\n")
      scale_silent_muts <- TRUE
      obs_exp_mut_count.df <- scale_exp_to_silent_muts(all_mutations.gr,
                                                       effect_filtered_mutations.gr,
                                                       test_regions.gr,
                                                       mut_effects.df,
                                                       mutfreqs.df,
                                                       config$silent_mut_flank_upstream,
                                                       config$silent_mut_flank_downstream)
      scaling_string <-  paste0("silent_bkg_mutrate_flank_up_", config$silent_mut_flank_upstream, "_down_", config$silent_mut_flank_downstream, "_")
      test_region_list <- test_region_list[test_region_list %in% filter(obs_exp_mut_count.df, silent_muts >= config$min_silent_mutations)$test_region]
    } 
  } else {
    cat("Incorrect seq_type in config. Quitting\n")
    quit()
  }
  # Set a scaling string for the output filename if no scaling was used.
  if(!(scale_with_regression || scale_silent_muts)){
    cat("No scaling\n")
    scaling_string <- "unadjusted_bkg_mutrate_"
  }
  
  
  ################################################
  # Filtering the test region list
  ################################################
  
  # Filter test regions to only include those with enough mutations
  test_region_list <- filter_test_regions_by_mut_count(test_region_list, effect_filtered_mutations.gr, config$min_mutations)
  
  if(config$remove_overlapping_test_regions){
    test_region_list <- test_region_list %>% filter_overlapping_test_regions(test_regions.gr)
  }
  
  if(config$plot_region_list_only){
    test_region_list <- test_region_list[test_region_list %in% config$plot_region_list]
  }
  
  # Make mut_effects.df a bit smaller by removing filtered test regions (after all filtering)
  mut_effects.df <- mut_effects.df %>% filter(test_region %in% test_region_list)
  
  # Split mut_effects.df into list of tibbles, so that we can pass them instead of the full mut_effects.df to foreach workers, to reduce memory usage.
  mut_effects.df <- split.data.frame(mut_effects.df, mut_effects.df$test_region)
  
  if(length(test_region_list) == 0){
    cat("No regions left for testing after filtering.\n")
    cat("Writing empty file...\n")
    
    empty.df <- data.frame(
        region = character(0),
        mutated_tumours = integer(0),
        exp_mutated_tumours = integer(0),
        p.cohort_skew = numeric(0),
        fdr.cohort_skew = numeric(0)
    )

    empty.df %>%
      write_tsv(paste0(config$out_dir,
                       config$out_base_name, ".", "seismic.tsv"))
    
    cat("Quitting\n")
    quit()
  }
  
  ################################################
  # Analysis loop
  ################################################

  # regular_output_mode FALSE only used for development, not in normal use.
  regular_output_mode <- TRUE
  # Remove unnecessary objects to decrease memory usage
  rm(all_mutations.gr, annotated_mutations.gr)
  gc()
  
  
  
  
  
  
  
  

  cat("Testing regions\n")
  pb <- txtProgressBar(max = length(test_region_list), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  cl <- makeCluster(no_of_cores, type="SOCK")
  clusterCall(cl, worker.init)
  registerDoSNOW(cl)
  
  ranks.df <- foreach(test_region_mut_effects = mut_effects.df,
                      .combine = 'bind_rows',
                      .errorhandling = 'remove',
                      .inorder = FALSE,
                      .multicombine = TRUE,
                      .options.snow = opts,
                      .packages = c('fitdistrplus', 'tidyverse', 'plyranges', 'data.table', 'stringi', 'inline', 'Rcpp', 'cowplot'),
                      .export = ls(.GlobalEnv)) %dopar% {

    current_test_region <- unique(test_region_mut_effects$test_region)
    
    
    # tryCatch the whole loop so that we can log any regions that encountered errors, instead of failing silently.
    tryCatch({ 
      
      plot_this_test_region <- current_test_region %in% config$plot_region_list
      
      
      test_region_mut_effects.gr <- test_region_mut_effects %>% as_granges()
      
      # Get trinuc info in test region.
      # For each base, get the trinuc, and include the effect of the 3 possible mutations
      test_region_mut_effects_tidy.df <- test_region_mut_effects.gr %>%
        as_tibble() %>%
        select(trinuc = ref_trinuc, C, A, G, T, base_no) %>%
        pivot_longer(-c(trinuc, base_no), names_to = "varnuc", values_to = "effect") %>% 
        filter(effect != "u") # 'u' refers to unchanged, e.g., A>A "mutation"
      
      
      
      # Negligible performance difference from defining the function inside the loop, and easier to read.
      calc_mutprobs_per_base <- function(test_region_mut_effects_tidy.df, combined_mutfreqs.df, effects_to_keep = NULL){
        # Further filter to effects included in config$effects_to_keep
        # Most common example: keep only non-silent mutations in a gene
        test_region_mut_effects_tidy_effect_filtered_combined_varnucs.df <- test_region_mut_effects_tidy.df %>% 
          { if(is.null(effects_to_keep)) . else filter(., effect %in% effects_to_keep) } %>% 
          group_by(trinuc, base_no) %>%
          arrange(varnuc) %>%
          summarise(combined_muttype = paste0(varnuc, collapse = ''), .groups = 'drop')
        
        # Get probability of mutation at each base in the test region
        mutprobs_per_base.df <- test_region_mut_effects_tidy_effect_filtered_combined_varnucs.df %>%
          inner_join(combined_mutfreqs.df, by = c('trinuc', 'combined_muttype')) %>% 
          select(sampleID, pmut = mutfreq)
        
        mutprobs_per_base.df
      }
      
      mutprobs_per_base.df <- calc_mutprobs_per_base(test_region_mut_effects_tidy.df, combined_mutfreqs.df, config$effects_to_keep)
      
  
      
      # Filter mutations to only look at those overlapping the test regions (e.g., the CDSs in a gene)
      mutations.df <- effect_filtered_mutations.gr %>% 
        filter_by_overlaps(
          test_region_mut_effects.gr %>% GenomicRanges::reduce() # Reduce probably unnecessary to perform manually.
          ) %>% 
        as_tibble()
      
      
  
      # Take each mutation in a each gene and make a tibble with that info as binary "mutated" status
      donor_mutated.df <- mutations.df %>% 
        group_by(sampleID) %>% 
        summarise(mutated = sign(dplyr::n()), .groups = 'drop')
      
      
      
      
      # Set scaling factor. For the cohort_dist results, i.e. when not taking recurrence into account, scaling here is irrelevant.
      if(scale_with_regression){
        # Have to redo the probability of mutation at each base to include all mutations (not filtering e.g. silent)
        # when calculating the expected number of mutations, as the observed mutation rate in reglist includes all mutations.
        no_of_exp_mutations <- calc_mutprobs_per_base(test_region_mut_effects_tidy.df, combined_mutfreqs.df) %>%
          pull(pmut) %>%
          sum()
        
        test_region_mutrate <- reglist$test_region_mutrate %>% filter(test_region == current_test_region)
        new_no_of_exp_mutations <- predict(reglist$lm, test_region_mutrate) * test_region_mutrate$length
        scaling_factor <- new_no_of_exp_mutations / no_of_exp_mutations
      } else if(scale_silent_muts){
        scaling_factor <- filter(obs_exp_mut_count.df, test_region == current_test_region)$obs_exp_silent_ratio
      } else {
        scaling_factor <- 1
      }
      # Mustn't set mutation probabilities negative, so abort loop iteration if scaling_factor<0
      if(scaling_factor < 0){
        # I would like a next statement here, but that doesn't seem to work with foreach. Therefore returning empty row.
        return(tibble(test=character(),
                      test_value=double(),
                      rank=double(),
                      region=character(),
                      obs_mutated_count=double(),
                      exp_mutated_count=double(),
                      scaling_factor=double()))
      }
      
      # Calculate the likelihood of any mutations (with the right effect) in the test region per patient,
      # scaling probabilities if so configured. Also add mutated status.
      mutprobs_test_region.df <- mutprobs_per_base.df %>%
        group_by(sampleID) %>% 
        summarise(exp_mut_count = sum(pmut * scaling_factor), pmut = 1 - prod(1 - pmut*scaling_factor), .groups = 'drop') %>% # 
        left_join(donor_mutated.df, by = "sampleID") %>% 
        replace_na(list(mutated = 0)) %>% 
        mutate(p_no_mut = 1 - pmut) %>%
        arrange(pmut)
      
      obs_no_of_mutated_patients <- mutprobs_test_region.df$mutated %>% sum()
      
      output.df <- tibble()
      if(config$test_recurrence_alone || config$test_combined_model){
        
        # We can use the same simulations for testing with only recurrence, and with the combined model
        sim_test_values <- sapply(1:config$no_of_simulations, function(x) sim_recurrence_and_combined_test_values(no_of_patients, mutprobs_test_region.df$pmut))
        
        
        if(config$test_recurrence_alone){
          
          # This is needed for the gamma distribution fit, since 0 is not allowed.
          # Can't be something really small like 10^-6 as that skews the distribution a lot, but 1 seems good.
          extra_mut <- 1
          
          fit.gamma <- fitdist(sim_test_values[1,] + extra_mut, distr = "gamma", method = "mle", lower = 0)
          
          output.df <- tibble(scaling_factor = scaling_factor,
                              expected_mutated_after_cohortdist_scaling = NA,
                              test_value = obs_no_of_mutated_patients,
                              rank = NA,
                              p = pgamma(obs_no_of_mutated_patients + extra_mut, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"], lower.tail = F),
                              test = "recurrence") %>% 
            bind_rows(output.df)
          
          
          
          
          
          
        }
        
        if(config$test_combined_model){
  
          
          test_value <- cohort_outcome_prob(mutprobs_test_region.df$p_no_mut, mutprobs_test_region.df$mutated)
          rank <- test_value %>% get_obs_rank(sim_test_values[2,])
          
          output.df <- tibble(scaling_factor = scaling_factor,
                              expected_mutated_after_cohortdist_scaling = NA,
                              test_value = test_value,
                              rank = rank,
                              p = rank / (config$no_of_simulations + 1),
                              test = "combined"
                              ) %>% 
            bind_rows(output.df)
          
        }
        
        
      }
      
      
      
      
      
      
      
      
      
      # This is the method that tests only for which patients have mutations, while ignoring recurrence
      if(config$test_cohort_skew_alone){
  
        # To help optimize() find the scaling factor that gets exp_mutated = obs_mutated, we start by estimating what the factor should be, ignoring the non-linearity of the relationship between p(mutated) and exp_no_of_mutations
        orig_exp_no_of_mutated_patients <- mutprobs_per_base.df %>%
          group_by(sampleID) %>% 
          summarise(pmut = 1 - prod(1 - pmut)) %>%
          pull(pmut) %>%
          sum()
        naive_scaling_factor <- obs_no_of_mutated_patients / orig_exp_no_of_mutated_patients
  
        # Try to find a scaling value where the expected number of MUTATED tumours is the same as what's observed.
        scaling_factor_cohortdist <- optimize(function(x) calc_exp_mutated_with_scaling(x, mutprobs_per_base.df, obs_no_of_mutated_patients), interval = c(0, naive_scaling_factor * 10), tol = 0.01)$minimum
        
        
        
        mutprobs_test_region_cohortdist.df <- mutprobs_per_base.df %>%
          group_by(sampleID) %>% 
          summarise(exp_mut_count = sum(pmut * scaling_factor_cohortdist), pmut = 1 - prod(1 - pmut*scaling_factor_cohortdist), .groups = 'drop') %>% # 
          left_join(donor_mutated.df, by = "sampleID") %>% 
          replace_na(list(mutated = 0)) %>% 
          mutate(p_no_mut = 1 - pmut) %>%
          arrange(pmut)
        
        
        expected_after_cohortdist_scaling <- sum(mutprobs_test_region_cohortdist.df$pmut)
        
        cohort_dist_scaling_ok <- TRUE
        if(abs(obs_no_of_mutated_patients - expected_after_cohortdist_scaling) > 1){
          cat("Cohort dist optimize scaling yields difference between observed and expected mutated patients of ",
              obs_no_of_mutated_patients - expected_after_cohortdist_scaling,
              " for ",
              current_test_region,
              ". Investigate whether the optimize function is finding the right scaling factor. \n")
          cohort_dist_scaling_ok <- FALSE
        }
        
  
        # With 0 mutations, the result will always be the same, so no point in running simulations. Faster to just output the answer (p -> 0.5 as no_of_sims -> Inf)
        # Normally, the tool should not be analysing genes with 0 mutations, but it could happen during troubleshooting, for instance.
        if(obs_no_of_mutated_patients == 0){
          output.df <- tibble(
            test_value = cohort_outcome_prob(mutprobs_test_region_cohortdist.df$p_no_mut, mutprobs_test_region_cohortdist.df$mutated),
            scaling_factor = scaling_factor_cohortdist,
            rank = config$no_of_simulations / 2 + 1,
            test = "cohort_dist") %>%
            bind_rows(output.df)
        } else {
          
          if(cohort_dist_scaling_ok){
            test_value <- cohort_outcome_prob(mutprobs_test_region_cohortdist.df$p_no_mut, mutprobs_test_region_cohortdist.df$mutated)
            
            # This doesn't work well for very low mutated counts, particularly just 1 (since that's not unimodal at all). Even 2 seems fine though.
            invisible(capture.output(
            fit.gamma <- fitdist(-sapply(1:config$no_of_simulations,
                                        function(x) cohort_outcome_prob(mutprobs_test_region_cohortdist.df$p_no_mut,
                                                                        sim_fixed_mutated_count2(no_of_patients, obs_no_of_mutated_patients, mutprobs_test_region_cohortdist.df$exp_mut_count))),
                                 distr = "gamma",
                                 method = "mle")
            ))
         
            p_val <- pgamma(-test_value, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"], lower.tail = F)
            output.df <- tibble(
              test_value = test_value,
              scaling_factor = scaling_factor_cohortdist,
              expected_mutated_after_cohortdist_scaling = expected_after_cohortdist_scaling,
              rank = NA,
              p = p_val,
              test = "cohort_dist") %>%
              bind_rows(output.df)
            
            # Plotting
            # If there is an error in the plotting, catch it so that only the plot and not the whole region result is affected.
            e <- tryCatch({
              if(p_val < config$plotting_max_p){plot_this_test_region <- TRUE}
              if(plot_this_test_region){
                p.cohort_dist_fit <- plot_cohort_dist_fit(fit.gamma, test_value, p_val)
                p.cohort_dist_cumsum <- plot_cohort_dist_cumsum(mutprobs_test_region_cohortdist.df, no_of_patients, obs_no_of_mutated_patients, config$effects_to_keep, config$no_of_simulations)
              }
            }, error = function(e) return())
            
          } else {
            # Output NA row if there was a problem getting the cohort dist scaling right.
            output.df <- tibble(
              test_value = NA,
              scaling_factor = scaling_factor_cohortdist,
              expected_mutated_after_cohortdist_scaling = expected_after_cohortdist_scaling,
              rank = NA,
              p = NA,
              test = "cohort_dist") %>%
              bind_rows(output.df)
          }
        }
        
      }
  
      # If no tests were performed, output dummy row so that we can get obs/exp values etc without running simulations
      if(nrow(output.df) == 0){
        output.df <- output.df %>% bind_rows(tibble(test = NA, test_value = NA, rank = NA))
      }
      
      
  
      # Plotting
      if(plot_this_test_region){
        # If there is an error in the plotting, catch it so that only the plot and not the whole region result is affected.
        e <- tryCatch({
          if(!exists('p.cohort_dist_fit') | !exists('p.cohort_dist_cumsum')){
            p.cohort_dist_fit <- NULL
            p.cohort_dist_cumsum <- NULL
          }
          p.combined <- plot_combined(
                        current_test_region,
                        config$cancer_type,
                        test_regions.gr,
                        mutations_with_effects_for_plotting.df,
                        config$reference_genome,
                        config$annotate_cds_effects,
                        no_of_patients,
                        p.cohort_dist_fit,
                        p.cohort_dist_cumsum,
                        highlight_position_min_muts = 2)
          pdf(file = paste0(fig_dir,
                            config$out_base_name, "_",
                            sig_string,
                            config$cancer_type,
                            "_", current_test_region, ".pdf"),
              width = 15,
              height = 15
              )
          print(p.combined)
          dev.off()
        }, error = function(e) return())
      }
      
      
      
      if(regular_output_mode){
        output.df %>%
          mutate(region = current_test_region,
                 obs_mutated_count = obs_no_of_mutated_patients,
                 exp_mutated_count = sum(mutprobs_test_region.df$pmut))
      } else {
        # This bit is used to look at patients' mutated status/mutations count, exp mut count and p_mutated for the cohort dist method
        mutprobs_test_region_cohortdist.df %>%
          left_join(mutations.df %>% count(sampleID, name = 'mutations')) %>%
          replace_na(list(mutations = 0)) %>%
          mutate(gene = current_test_region)
      }
      
      
    # If anything goes wrong inside the loop, return an ERROR row so that we know. 
    }, error = function(e) tibble(region = current_test_region, test = "ERROR"))  
    
  }
  # End of loop
  stopCluster(cl)
  cat("\n") # Just to get a newline after the progress bar.
  
  
  if("ERROR" %in% ranks.df$test){
    cat('Errors were encountered in the following regions (the relevant results will be missing):\n')
    error_regions <- ranks.df %>% 
      filter(test == 'ERROR') %>% 
      pull(region) %>% 
      unique()
    cat(paste0('  ', error_regions), sep = '\n')
  }
  
  ##################################
  # Output
  ##################################
  if(!config$plot_region_list_only){
    col_rename <- c(mutated_tumours = 'obs_mutated_count',
                    exp_mutated_tumours = 'exp_mutated_count',
                    p.cohort_skew = 'p.cohort_dist',
                    fdr.cohort_skew = 'fdr.cohort_dist',
                    p.skew_and_recurrence = 'p.combined',
                    fdr.skew_and_recurrence = 'fdr.combined')
    ranks.df %>%
      filter(test != 'ERROR') %>% 
      arrange(test, p) %>%
      group_by(test) %>% 
      mutate(fdr = p.adjust(p, method = 'fdr')) %>% 
      select(region, test, p, fdr, obs_mutated_count, exp_mutated_count) %>% 
      pivot_wider(names_from = test, values_from = c(p, fdr), names_sep = ".") %>% 
      dplyr::rename(any_of(col_rename)) %>% 
      write_tsv(paste0(config$out_dir,
                       config$out_base_name, "_",
                       scaling_string,
                       sig_string,
                       config$cancer_type,
                       "_min_", config$min_mutations, "_muts",
                       "_", config$no_of_simulations, "_sims.tsv"))
  }
    
  
  cat("Done!\n")
}






































################################################
# Read/write files
################################################


# Read SNVs from MAF file. If there is a cancer type column (name specified in cancer_colname), filter the mutations to only the ones matching the cancer types in config.
# If there is no cancer type column, make one and set it to the cancer type in config (just to be compatible with other parts of the script).
read_mutations_from_maf <- function(path, genome, cancer_type, patient_colname = 'Tumor_Sample_Barcode',  cancer_colname = 'cancer_colname_not_available'){
  patient_colname <- tolower(patient_colname)
  cancer_colname <- tolower(cancer_colname)
  bases <- c('C', 'A', 'G', 'T')
  pyrimidines <- c('C', 'T')
  df <- data.table::fread(path) %>% 
    as_tibble() %>% 
    rename_all(tolower) # PCAWG file doesn't follow the format conventions (e.g., End_position instead of End_Position), so let's just lc everything.
  
  if(!(patient_colname %in% colnames(df))){
    cat('patient column name not found in .maf file. Incorrectly spelled in the config file?\n')
    quit()
  }
  if(cancer_colname != 'cancer_colname_not_available' & !(cancer_colname %in% colnames(df))){
    cat('cancer column name not found in .maf file. Incorrectly spelled in the config file?\n')
    quit()
  }
  
  df %>% 
    { if( length(cancer_type) == 1 & !(cancer_colname %in% names(.)) ) mutate(., cancer = cancer_type) else if( cancer_colname %in% names(.) ) dplyr::rename(., cancer = (!!sym(cancer_colname)))} %>% # Cancer type column handling - add if missing.
    { if( length(cancer_type) == 1 & cancer_type[1] == 'pancancer' ) . else filter(., cancer %in% cancer_type) } %>% # Cancer type filtration
    filter(variant_type == 'SNP') %>% 
    mutate(varnuc = ifelse(reference_allele != tumor_seq_allele1, tumor_seq_allele1, tumor_seq_allele2)) %>% 
    dplyr::rename(seqnames = chromosome, start = start_position, end = end_position, old_strand = strand, refnuc = reference_allele, sampleID = (!!sym(patient_colname))) %>% 
    select(seqnames, start, end, old_strand, refnuc, varnuc, sampleID, cancer) %>% 
    mutate(seqnames = ifelse(!str_detect(seqnames, '^chr'), paste0('chr', seqnames), seqnames)) %>% 
    mutate(old_strand = ifelse(as.character(old_strand) == '1', '+', old_strand),
           old_strand = ifelse(as.character(old_strand) == '-1', '-', old_strand),
           strand = ifelse(refnuc %in% pyrimidines, old_strand, chartr('+-', '-+', old_strand))) %>% 
    mutate(refnuc = ifelse(strand == old_strand, refnuc, str_revcomp(refnuc)), 
           varnuc = ifelse(strand == old_strand, varnuc, str_revcomp(varnuc))) %>% 
    as_granges() %>% 
    mutate(trinuc = getSeq(genome, . + 1, as.character = TRUE)) %>% 
    select(cancer, sampleID, varnuc, trinuc)
}

# Read mutations from .tsv file
# filter for correct cancer type, and return GRanges with sampleID, varnuc, and trinuc
read_mutations_from_tsv <- function(mutations_path, genome, cancer_type){
  pyrimidines <- c('C', 'T')
  bases <- c('C', 'A', 'G', 'T')
  read_tsv(mutations_path, col_types = cols(refnuc = col_character(), varnuc = col_character())) %>% # Specifying col_type for ref/varnuc to avoid T being interpreted as TRUE
    { if(length(cancer_type) == 1 & cancer_type[1] == 'pancancer') . else filter(., cancer %in% cancer_type) } %>%
    dplyr::rename(old_strand = strand) %>% 
    filter(refnuc %in% bases, varnuc %in% bases) %>% 
    mutate(strand = ifelse(refnuc %in% pyrimidines, old_strand, chartr('+-', '-+', old_strand))) %>% 
    mutate(refnuc = ifelse(strand == old_strand, refnuc, str_revcomp(refnuc)), 
           varnuc = ifelse(strand == old_strand, varnuc, str_revcomp(varnuc))) %>% 
    as_granges() %>%
    mutate(trinuc = getSeq(genome, . + 1, as.character = TRUE)) %>% 
    select(cancer, sampleID, varnuc, trinuc)
}

# Used to check if the mutation file is our legacy format
is_legacy_mutation_format <- function(path, cancer_type){
  columns <- fread(path, nrows = 0) %>% colnames()
  all(c("seqnames", "start", "end", "strand", "cancer", "refnuc", "varnuc", "sampleID") %in% columns)
}

# Read mutations with read_mutations_from_tsv if the mutation file is in our legacy format. Otherwise, assume MAF.
read_mutations <- function(path, cancer_type, genome, patient_colname = 'Tumor_Sample_Barcode',  cancer_colname = 'cancer_colname_not_available'){
  if(is_legacy_mutation_format(path)){
    mutations <- read_mutations_from_tsv(path, genome, cancer_type)
  } else {
    mutations <- read_mutations_from_maf(path, genome, cancer_type, patient_colname, cancer_colname)
  }
  return(mutations)
}

# Read test regions, either with or with genome assembly filtering (hg19/hg38)
read_test_regions <- function(test_regions_path, assembly_arg){
  con <- file(test_regions_path, "r")
  colnames_in_file <- readLines(con, n=1) %>% str_split('\t') %>% unlist()
  close(con)
  
  assembly_both_formats <- get_chromosome_synonyms(assembly_arg)
  
  exp_colnames_no_assembly <- c('seqnames', 'start', 'end', 'region', 'strand')
  
  col_types <- list(seqnames = col_character(), start = col_integer(), end = col_integer(), strand = col_character(), region = col_character())
  if(setequal(colnames_in_file, exp_colnames_no_assembly)){
    cat("  Loading all regions in", test_regions_path, "\n")
    test_regions <- read_tsv(test_regions_path, col_types = do.call(cols, col_types))
  } else if(setequal(colnames_in_file, c(exp_colnames_no_assembly, 'assembly'))){
    cat("  Loading", assembly_arg, "regions in", test_regions_path, "\n")
    test_regions <- read_tsv(test_regions_path, col_types = do.call(cols, c(col_types, list(assembly = col_character())))) %>% 
      filter(assembly %in% assembly_both_formats) %>% 
      select(-assembly)
  } else {
    cat("  Column names not as expected in ", test_regions_path, "\n")
    quit()
  }
  
  test_regions %>% 
    dplyr::rename(test_region = region) %>% 
    as_granges() %>%
    arrange(test_region, seqnames, start)
}



################################################
# Mutation annotation and effect filtering
################################################

save_mut_effects <- function(mut_effects.df, test_region_path, rds_path, md5_path){
  cat("Saving annotated regions, to speed up next SEISMIC run with this region file\n")
  tryCatch({
    saveRDS(mut_effects.df, rds_path)
    system(paste('cd', dirname(test_region_path),'; md5sum', basename(test_region_path), basename(rds_path), '>', basename(md5_path)))
  }, 
  error = function(e) cat("  Couldn't save annotated regions - Annotation will have to be done again next run\n"))
}


get_mut_effects <- function(test_regions.gr, genome, assembly, test_region_path, annotate_cds_effects, no_of_cores){
  new_name <- paste0(test_region_path, '.annotated_regions', ifelse(annotate_cds_effects, '_with_mut_effects', ''), '_', assembly)
  rds_path <- paste0(new_name, '.rds')
  md5_path <- paste0(new_name, '.md5')
  if(all(file.exists(test_region_path, rds_path, md5_path))){
    tryCatch(command_status <- system(paste('cd', dirname(test_region_path), '; md5sum -c --status', basename(md5_path))),
             error = function(e) NULL)
    load_mut_effects <- command_status == 0
  } else {
    load_mut_effects <- FALSE
  }
  
  if(load_mut_effects){
    cat("Loading previously annotated regions\n")
    mut_effects.df <- readRDS(rds_path)
  } else {
    # If annotate_cds_effects is TRUE, annotate test_regions and mutations m/s/n, and filter by effects_to_keep
    # If annotate_cds_effects is FALSE, annotate test_regions and mutations na, and keep all mutations, regardless of what effects_to_keep is set to.
    if(annotate_cds_effects){
      cat("Annotating regions with mutation effects\n")
    } else {
      cat("Annotating regions without mutation effects\n")
    }
    codon_table.df <- make_codon_table()
    mut_effects.df <- do.call(bind_rows, mclapply(unique(test_regions.gr$test_region), function(x) annotate_mut_effects(test_regions.gr, x, genome, codon_table.df, annotate_cds_effects), mc.cores = no_of_cores))
    
    save_mut_effects(mut_effects.df, test_region_path, rds_path, md5_path)
  }
  
  return(mut_effects.df)
}

annotate_mut_effects <- function(cds.gr, gene_val, genome, codon_table.df, annotate_cds_effects){
  pyr <- c("C", "T")
  # Tile into 1 bp segments and get sequence. Enumerate codons and bases
  tiled_cds.df <- cds.gr %>%
    filter(test_region == gene_val) %>% 
    arrange(start) %>%
    tile_ranges(1) %>%
    select(-partition) %>%
    mutate(ref_trinuc = getSeq(genome, . + 1, as.character = TRUE),
           refnuc = str_sub(ref_trinuc, 2, 2)) %>%
    as_tibble() %>% 
    mutate(base_no = ifelse(strand == "+", row_number(), -(row_number() - dplyr::n() - 1)),
           codon_no = ceiling(base_no / 3))

  if(annotate_cds_effects){
    region_strand <- unique(tiled_cds.df$strand) %>% as.character()
    region_seq <- ifelse(region_strand == "+", paste0(tiled_cds.df$refnuc, collapse = ""), reverse(paste0(tiled_cds.df$refnuc, collapse = "")))
    
    # Annotate each possible mutation at each base with effect
    gene_codons.df <- tibble(codon = strsplit(region_seq, "(?<=.{3})", perl = TRUE)[[1]]) %>% 
      mutate(codon_no = row_number()) %>% 
      dplyr::slice(rep(1:dplyr::n(), each = 12)) %>% 
      mutate(codon_base_no = rep(1:3, dplyr::n()/12, each = 4)) %>% 
      mutate(mut_to = rep(c("A", "C", "G", "T"), dplyr::n()/4)) %>% 
      mutate(codon1 = substr(codon, 1, codon_base_no - 1),
             codon2 = substr(codon, codon_base_no + 1, 3)) %>% 
      mutate(mut_codon = paste0(codon1, mut_to, codon2)) %>% 
      select(-codon1, -codon2) %>% 
      left_join(codon_table.df, by = 'codon') %>% 
      left_join(codon_table.df %>% dplyr::rename(mut_codon = codon, mut_aa = aa), by = 'mut_codon') %>% 
      mutate(effect = case_when(codon == mut_codon ~ "u", aa == mut_aa ~ "s", mut_aa == "STOP" ~ "n", TRUE ~ "m"))
    # u - unmutated
    # s - synonymous
    # n - nonsense
    # m - missense
    
    # Sorting mutation effects and pivoting wider to be able to bind_cols to tiled_cds.df
    # Making pyrimidine-focused
    gr_sorted_effects.df <-
      gene_codons.df %>%
      mutate(refnuc = substr(codon, codon_base_no, codon_base_no)) %>%
      mutate(strand_switch = ifelse(refnuc %in% pyr, FALSE, TRUE)) %>%
      mutate(mut_to = ifelse(strand_switch, str_revcomp(mut_to), mut_to),
             refnuc = ifelse(strand_switch, str_revcomp(refnuc), refnuc)) %>%
      select(codon_no, codon_base_no, mut_to, effect, strand_switch) %>%
      pivot_wider(names_from = mut_to, values_from = effect) %>% 
      mutate(strand = region_strand) %>% 
      mutate(plus_strand_codon_order = ifelse(strand == "+", codon_no, -(codon_no - dplyr::n()/3 - 1)),
             plus_strand_codon_base_order = ifelse(strand == "+", codon_base_no, -(codon_base_no - 4))) %>% 
      arrange(plus_strand_codon_order, plus_strand_codon_base_order)
    
    # Adding mutation effects to tiled cds
    tiled_cds.df %>%
      ungroup() %>%
      bind_cols(gr_sorted_effects.df %>% select(C,A,G,T, strand_switch)) %>%
      mutate(strand = ifelse(strand_switch, chartr("+-", "-+", as.character(strand)), as.character(strand)),
             refnuc = ifelse(strand_switch, str_revcomp(refnuc), refnuc),
             ref_trinuc = ifelse(strand_switch, str_revcomp(ref_trinuc), ref_trinuc),
             test_region = gene_val) %>% 
      select(-width, -strand_switch)
  } else {
    tiled_cds.df %>% mutate(strand_switch = ifelse(refnuc %in% pyr, FALSE, TRUE),
                            strand = ifelse(strand_switch, chartr("+-", "-+", as.character(strand)), as.character(strand)),
                            refnuc = ifelse(strand_switch, str_revcomp(refnuc), refnuc),
                            ref_trinuc = ifelse(strand_switch, str_revcomp(ref_trinuc), ref_trinuc),
                            C = 'na',
                            A = 'na',
                            G = 'na',
                            T = 'na',
                            test_region = gene_val) %>%
      select(-width, -strand_switch) %>% 
      pivot_longer(c(C, A, G, T), values_to = 'effect', names_to = 'varnuc') %>%
      mutate(effect = ifelse(refnuc == varnuc, 'u', effect)) %>%
      pivot_wider(names_from = varnuc, values_from = effect)
  }
}

# Annotate mutations in mut_effects regions
annotate_mutations <- function(mutations.gr, mut_effects.df){
  mutations.gr %>%
    as_tibble() %>%
    mutate(seqnames = as.character(seqnames), strand = as.character(strand)) %>% # For joining. Might be better to make mut_effects columns factors, but I'm wary of getting the levels wrong
    inner_join(mut_effects.df, by = c("seqnames", "start", "end", "strand")) %>%
    pivot_longer(c("C","A","G","T"), names_to = "effect_varnuc", values_to = "effect") %>% 
    filter(effect_varnuc == varnuc) %>%
    select(-effect_varnuc) %>% 
    as_granges()
}

# Filter mutation effect by what to keep
filter_mutation_effect_keep <- function(mutations_with_effects.gr, effects){
  filter(mutations_with_effects.gr, effect %in% effects)
}

# Filter mutation effect by what to remove
filter_mutation_effect_remove <- function(mutations_with_effects.gr, effects){
  filter(mutations_with_effects.gr, !(effect %in% effects))
}

make_codon_table <- function(){
  tribble(~aa, ~codon,
          'G', 'GGG',
          'G', 'GGA',
          'G', 'GGT',
          'G', 'GGC',
          'E', 'GAG',
          'E', 'GAA',
          'D', 'GAT',
          'D', 'GAC',
          'V', 'GTG',
          'V', 'GTA',
          'V', 'GTT',
          'V', 'GTC',
          'A', 'GCG',
          'A', 'GCA',
          'A', 'GCT',
          'A', 'GCC',
          'R', 'AGG',
          'R', 'AGA',
          'S', 'AGT',
          'S', 'AGC',
          'K', 'AAG',
          'K', 'AAA',
          'N', 'AAT',
          'N', 'AAC',
          'M', 'ATG',
          'I', 'ATA',
          'I', 'ATT',
          'I', 'ATC',
          'T', 'ACG',
          'T', 'ACA',
          'T', 'ACT',
          'T', 'ACC',
          'W', 'TGG',
          'STOP', 'TGA',
          'C', 'TGT',
          'C', 'TGC',
          'STOP', 'TAG',
          'STOP', 'TAA',
          'Y', 'TAT',
          'Y', 'TAC',
          'L', 'TTG',
          'L', 'TTA',
          'F', 'TTT',
          'F', 'TTC',
          'S', 'TCG',
          'S', 'TCA',
          'S', 'TCT',
          'S', 'TCC',
          'R', 'CGG',
          'R', 'CGA',
          'R', 'CGT',
          'R', 'CGC',
          'Q', 'CAG',
          'Q', 'CAA',
          'H', 'CAT',
          'H', 'CAC',
          'L', 'CTG',
          'L', 'CTA',
          'L', 'CTT',
          'L', 'CTC',
          'P', 'CCG',
          'P', 'CCA',
          'P', 'CCT',
          'P', 'CCC')
}

################################################
# Mutation frequencies
################################################

# Function to calculate mutation frequencies per trinuc and varnuc in each patient.
# For WXS data, estimate patient-specific mutation frequencies by scaling cohort-level mutation frequencies to mutational burden in patients
calculate_mutfreqs <- function(mutations.gr, sig_type, seq_type, genome, assembly_arg, wxs_mutfreq_regions){
  # Get trinucleotide counts in exons (WXS) or chr1-22 (WGS). Limit mutations to the same regions for the mutfreq calculation.
  if(seq_type == 'WGS'){
    gr <- get_genome_gr(genome, paste0('chr', 1:22))
  } else if(seq_type == 'WXS'){
    # If WXS data, filter the exons data file to get the right assembly exon definitions
    assembly_both_formats <- get_chromosome_synonyms(assembly_arg)
    gr <- read_tsv(wxs_mutfreq_regions,
                   col_types = cols(
                     seqnames = col_character(),
                     start = col_integer(),
                     end = col_integer(),
                     assembly = col_character()
                   )) %>%
      filter(assembly %in% assembly_both_formats) %>%
      select(-strand) %>% # Get rid of strand before granges conversion so that reduce doesn't leave overlaps on different strands.
      as_granges() %>% 
      GenomicRanges::reduce()
  } else{
    cat("Incorrect seq_type in config. Quitting\n")
    quit()
  }
  trinuc_bg_counts.df <- get_genome_trinuc_bgfreqs(genome, gr)
  mutations.gr <- mutations.gr %>% filter_by_overlaps(gr)
  
  # All combinations of patient, trinuc (pyr-centric) and varnuc (in order to add missing rows to mut_burden/mut_freqs)
  DNA_bases <- c("A", "C", "G", "T")
  patients <- mutations.gr$sampleID %>% unique()
  all_combos <- expand_grid(varnuc=DNA_bases, trinuc=get_pyr_trinucs(), sampleID=patients) %>%
    filter(substr(trinuc, 2, 2) != varnuc)
  
  # Separate mutfreq calculations for wgs and wxs, since wxs doesn't have enough mutations to do it the same way as wgs. Can use the wxs method for wgs as well (use_cohort_sig)
  if(sig_type == "cohort"){
    # Calculate mutation counts per patient 
    mut_burden.df <- mutations.gr %>%
      as_tibble() %>%
      count(cancer, sampleID) %>%
      arrange(desc(n))
    
    muts_per_cancer_type.df <- mutations.gr %>% 
      as_tibble() %>% 
      count(cancer, name = 'muts_per_cancer_type')
    
    # Calculate mutfreq for all patients together, and normalise by mutation count to get mutfreq_per_mut
    cohort_mutfreq.df <- mutations.gr %>%
      as_tibble() %>%
      count(cancer, varnuc, trinuc, name = "burden") %>%
      left_join(trinuc_bg_counts.df, by = "trinuc") %>%
      left_join(muts_per_cancer_type.df, by = 'cancer') %>% 
      mutate(mutfreq = burden / count) %>% 
      mutate(mutfreq_per_mut = mutfreq / muts_per_cancer_type)
    
    # Combine cohort mutfreq and patient mutation burdens to estimate patient mutfreqs.
    mutfreqs.df <- all_combos %>% 
      left_join(mut_burden.df, by = "sampleID") %>% 
      left_join(cohort_mutfreq.df %>% select(cancer, trinuc, varnuc, mutfreq_per_mut), by = c("cancer", "varnuc", "trinuc")) %>% 
      mutate(mutfreq = mutfreq_per_mut * n) %>% 
      select(-mutfreq_per_mut, -n)
    
  } else if(sig_type == "patient"){
    # Calculate mutation burden per trinuc and varnuc, for each patient
    mut_burden.df <- mutations.gr %>%
      as_tibble() %>%
      count(cancer, sampleID, varnuc, trinuc, name = "burden") %>%
      complete(all_combos) %>%
      replace_na(list(burden = 0))
    
    # Combine with trinuce bg counts and calculate mutfreqs
    mutfreqs.df <- mut_burden.df %>% 
      left_join(trinuc_bg_counts.df, by = "trinuc") %>% 
      mutate(mutfreq = burden / count) %>%
      select(-burden, -count)
    
  } else if(sig_type == "flat"){
    mut_burden.df <- mutations.gr %>% 
      as_tibble() %>% 
      count(cancer, sampleID) %>% 
      arrange(desc(n))
    
    total_trinucs <- trinuc_bg_counts.df %>% filter(substr(trinuc, 2, 2) %in% c("C", "T")) %>% pull(count) %>% sum()
    
    mutfreqs.df <- mut_burden.df %>%
      mutate(mutfreq = n / total_trinucs) %>% 
      left_join(all_combos, by = "sampleID") %>% 
      mutate(mutfreq = mutfreq / 3) %>% # Equal distribution between the three possible mutations at a base
      select(-n)
  } else {
    cat("wrong sig_type\n")
    quit()
  }
  
  # Make sure mutation types that aren't represented are 0 instead of NA when returning
  mutfreqs.df %>% 
    replace_na(list(mutfreq = 0))
}

# Expand mutfreqs df to show the sum of each combination of varnucs for each trinuc.
# Basically useful so that we can calculate the expected number of mutations for each patient at any site, and only do it once instead of once per gene.
# Example, what's the mutfreq for patient n at a TCC>G/A (TCC>T silent in this scenario)? Useful.
combine_mutfreq_varnucs <- function(mutfreqs.df){
  combined_mutfreqs.df <- tibble()
  for(trinuc_iter in get_pyr_trinucs()){
    for(combined_muttype_iter in c("A", "C", "G", "T", "AC", "AG", "AT", "CG", "CT", "GT", "ACG", "ACT", "AGT", "CGT")){# Viktigt att de är i bokstavsordning (AGT istället för ATG)
      separate_muttypes <- combined_muttype_iter %>% str_split('') %>% unlist()
      if( !(str_sub(trinuc_iter, 2, 2) %in% separate_muttypes)){
        combined_mutfreqs.df <- mutfreqs.df %>%
          filter(trinuc == trinuc_iter, varnuc %in% separate_muttypes) %>% 
          group_by(sampleID, cancer, trinuc) %>% 
          summarise(mutfreq = sum(mutfreq), .groups = 'drop') %>% 
          mutate(combined_muttype = combined_muttype_iter) %>% 
          bind_rows(combined_mutfreqs.df)
      }
    }
  }
  return(combined_mutfreqs.df)
}





################################################
# Scaling mutation expectations - WXS
################################################

initiate_regression_list <- function(test_regions.gr, all_mutations.gr){
  reglist <- list()
  reglist$variables <- character()
  reglist$scale_with_regression_string <- character()
  reglist$test_region_mutrate <- calculate_test_region_mutrate(test_regions.gr, all_mutations.gr)
  return(reglist)
}

calculate_test_region_mutrate <- function(test_regions.gr, all_mutations.gr){
  # Get test region lengths
  test_region_lengths.df <- test_regions.gr %>%
    as_tibble() %>%
    group_by(test_region) %>%
    summarise(length = sum(width), .groups = 'drop')
  
  # Calculate mutation rates in test regions
  hits <- findOverlaps(test_regions.gr, all_mutations.gr, ignore.strand = TRUE)
  mutrate.df <- tibble(test_region = test_regions.gr[hits@from]$test_region) %>% 
    count(test_region, name = "mutations") %>% 
    inner_join(test_region_lengths.df, by = "test_region") %>% 
    mutate(mutrate = mutations / length)
  
  return(mutrate.df)
}

# Add replication timing data and variable to reglist
scale_exp_to_repl_reg <- function(reglist, repl_timing_path, test_regions.gr){
  # Load repl timing data
  repl_timing.gr <- import.bed(repl_timing_path) %>% 
    filter(score != 0) # 0 means undefined, I think. Have to get rid of it.
  
  # Assign S50 values to genes
  hits <- findOverlaps(test_regions.gr, repl_timing.gr, ignore.strand = TRUE)
  test_region_timing.df <- tibble(test_region = test_regions.gr[hits@from]$test_region, S50 = repl_timing.gr[hits@to]$score) %>% 
    group_by(test_region) %>% 
    summarise(S50 = mean(S50), .groups = 'drop')
  
  reglist$test_region_mutrate <- reglist$test_region_mutrate %>% inner_join(test_region_timing.df, by = "test_region")
  reglist$variables <- c(reglist$variables, "S50")
  reglist$scale_with_regression_string <- c(reglist$scale_with_regression_string, "repl_timing_")
  
  return(reglist)
}

# Add expression data and variable to reglist
scale_exp_to_expr_reg <- function(reglist, expression_path, cancer_type){
  expression.df <- read_tsv(expression_path, col_types = cols(gene = col_character(), expression = col_double(), cancer = col_character())) %>%
    filter(cancer == cancer_type) %>% 
    filter(expression != 0) %>% # Necessary if we use log(expr). Removes a lot of genes, but not many that would otherwise be included, I think.
    mutate(log_expr = log(expression))
  
  if(nrow(expression.df) == 0){
    return(reglist)
  }
  
  reglist$test_region_mutrate <- reglist$test_region_mutrate %>% inner_join(expression.df %>% dplyr::rename(test_region = gene), by = "test_region") # This renaming might not work well if we're analysing something other than genes - we still want to do this operation by gene.
  reglist$variables <- c(reglist$variables, "log_expr")
  reglist$scale_with_regression_string <- c(reglist$scale_with_regression_string, "expression_")
  
  return(reglist)
}

# Make regression with data and variables in reglist
make_scaling_regression <- function(reglist, min_test_region_length_for_regression){
  f <- as.formula(
    paste("mutrate", 
          paste(reglist$variables, collapse = " + "), 
          sep = " ~ "))
  
  reglist$lm <- reglist$test_region_mutrate %>% 
    filter(length > min_test_region_length_for_regression) %>% 
    lm(f, data = .)
  reglist$scale_with_regression_string <- c(reglist$scale_with_regression_string,
                                            "bkg_mutrate_min_",
                                            min_test_region_length_for_regression,
                                            "_bp_genes_reg_") %>% 
    paste0(collapse="")
  
  return(reglist)
}





################################################
# Scaling mutation expectations - WGS
################################################

# Count silent mut positions in test_regions.gr, and around them.
# Note that only the effect "s" counts as silent, so if we're analysing e.g. TFBSs, then we can avoid counting
# the mutations in them by annotating them with anything else, like "na". Then, only the mutations around them will be used.
count_silent_mut_positions <- function(mut_effects.df, test_regions.gr, flank_upstream, flank_downstream){
  # Make GRanges with one range per test region, encompassing all ranges for each test region (e.g., make one range with the transcript based on cds regions)
  # Then extend upstream and downstream by desired amounts.
  # Requires that all ranges in a test_region are on the same chromosome and strand
  test_region_span_with_flanks.gr <- test_regions.gr %>% 
    get_gr_region_range("test_region") %>% 
    extend_gr_flanks(flank_upstream, flank_downstream)

  # Remove parts of test_region_span_with_flanks.gr that overlap with test_regions.gr TEST_REGION-WISE.
  # This means that transcribed regions (minus cds regions) are allowed to belong to multiple genes, for the purposes of counting synonymous mutations, which are used to estimate gene-level background mutation rate.
  # Since I can't group_by gene and use intersect on each group, I combined chrom and gene into artificial seqlevels, pretending each gene is on a different chromosome.
  # This hack means I can use intersect gene-wise, which is WAY faster than e.g. splitting up the txn.gr based on gene and using intersect once per gene.
  silent_regions.gr <- setdiff(test_region_span_with_flanks.gr %>% make_fake_paste_column_chrom("test_region", sep = "¨"),
                               test_regions.gr %>% make_fake_paste_column_chrom("test_region", sep = "¨")) %>%
    as_tibble() %>%
    separate(seqnames, into=c("chrom", "test_region"), sep = "¨") %>%
    select(seqnames = chrom, start, end, strand, test_region) %>%
    as_granges()
  
  # For each gene, count the number of trinucleotides in introns + gene flanks.
  # Use that count for each possible varnuc, as every SNV in an intron counts as silent
  intron_trinuc_count.df <- getSeq(genome, silent_regions.gr) %>% 
    trinucleotideFrequency() %>% 
    as_tibble() %>% 
    mutate(test_region = silent_regions.gr$test_region) %>%
    pivot_longer(-test_region, names_to="trinuc", values_to="count") %>% 
    group_by(test_region, trinuc) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    mutate(trinuc = ifelse(substr(trinuc, 2, 2) %in% c("C", "T"), trinuc, str_revcomp(trinuc))) %>% 
    group_by(test_region, trinuc) %>% 
    summarise(count = sum(count), .groups = "drop") %>% 
    mutate(A = count, T = count, C = count, G = count) %>% 
    select(-count) %>% 
    pivot_longer(c(A, T, C, G), names_to = "varnuc", values_to = "count") %>% 
    filter(substr(trinuc, 2, 2) != varnuc)
    
  # For each rest_region (gene), count the number of silent mutation positions in each combination of trinuc and varnuc
  # Current data.table solution much faster than previous dbplyr version.
  exon_silent_count.df <- mut_effects.df %>% 
    dplyr::rename(trinuc=ref_trinuc) %>% 
    setDT() %>% 
    .[,list(T=ifelse(T == "s", 1, 0),
            C=ifelse(C == "s", 1, 0),
            G=ifelse(G == "s", 1, 0),
            A=ifelse(A == "s", 1, 0),
            test_region,
            trinuc)] %>% 
    .[, list(T = sum(T), C = sum(C), G = sum(G), A = sum(A)), .(test_region, trinuc)] %>% 
    as_tibble() %>% 
    pivot_longer(c(C, A, G, T), names_to = "varnuc", values_to="count") %>% 
    filter(substr(trinuc, 2, 2) != varnuc)
    
  all_combos <- expand_grid(varnuc=c("A", "C", "G", "T"),
                            trinuc=get_pyr_trinucs(),
                            test_region = unique(test_regions.gr$test_region)) %>%
    filter(substr(trinuc, 2, 2) != varnuc) %>% 
    mutate(count = 0)
  
  silent_mut_pos.df <- exon_silent_count.df %>% 
    bind_rows(intron_trinuc_count.df) %>% 
    bind_rows(all_combos) %>% 
    group_by(test_region, trinuc, varnuc) %>% 
    summarise(count = sum(count), .groups = 'drop')
  
  return(silent_mut_pos.df)
}


scale_exp_to_silent_muts <- function(all_mutations.gr, effect_filtered_mutations.gr, test_regions.gr, mut_effects.df, mutfreqs.df, flank_upstream, flank_downstream){
  silent_mut_pos.df <- count_silent_mut_positions(mut_effects.df, test_regions.gr, flank_upstream, flank_downstream)
  
  # calculate the expected number of silent mutations per gene for all genes we've decided to analyse .
  expected_silent_muts.df <- silent_mut_pos.df %>%
    # filter(test_region %in% test_region_list) %>% 
    dplyr::rename(bg_count = count) %>%
    inner_join(mutfreqs.df, by = c("trinuc", "varnuc")) %>%
    mutate(exp_silent_muts = mutfreq * bg_count) %>% 
    setDT() %>% 
    .[, list( exp_silent_muts = sum(exp_silent_muts)), .(test_region)] %>% # This summarise could crash with too big datasets using tidyverse
    as_tibble()
  
  
  
  # Get mutations that are in effect_filtered_mutations. OBS! Not necessarily all non-silent. Depends on what's specified in effect_to_keep.
  # But anything other than that effect will be regarded as silent. Might need to change this at some point, but I only see myself using it
  # to either count non-silent mutations as non-silent, or all mutations as non-silent, in case we don't want to include any mutations in
  # the test_region in the silent mutation count (e.g. in a TFBS)
  nonsilent_mut_count.df <- effect_filtered_mutations.gr %>%
    get_gr_region_range_count("test_region", "nonsilent_muts")
  
  test_region_span_with_flanks.gr <- test_regions.gr %>% 
    get_gr_region_range("test_region") %>% 
    extend_gr_flanks(flank_upstream, flank_downstream)
  
  hits <- findOverlaps(all_mutations.gr, test_region_span_with_flanks.gr, ignore.strand = TRUE)
  mut_count_per_test_region.df <- tibble(test_region = test_region_span_with_flanks.gr$test_region[hits@to]) %>% 
    count(test_region, name = 'all_mutations') %>% 
    left_join(nonsilent_mut_count.df, by = "test_region") %>% 
    replace_na(list(nonsilent_muts = 0)) %>% 
    mutate(silent_muts = all_mutations - nonsilent_muts)
  
  # Calculate the ratio of observed and expected silent mutations, used to scale the expected number of mutations
  obs_exp_mut_count.df <- mut_count_per_test_region.df %>% 
    left_join(expected_silent_muts.df, by = "test_region") %>% 
    mutate(obs_exp_silent_ratio = silent_muts / exp_silent_muts)
  
  return(obs_exp_mut_count.df)
}




################################################
# Scaling mutation expectations - Cohort distribution
################################################

calc_exp_mutated_with_scaling <- function(scaling, mutprobs_per_base.df, observed_mutated_count){
  new_exp <- mutprobs_per_base.df %>%
    group_by(sampleID) %>% 
    summarise(pmut = 1 - prod(1 - pmut*scaling), .groups = 'drop') %>% 
    pull(pmut) %>% 
    sum()
  abs(log2(new_exp / observed_mutated_count))
}





################################################
# Filter test regions
################################################

# Filter test regions (e.g. genes) by how many mutations of a certain effect (e.g. non-silent - i.e. n and m) they have.
# Return list of regions that have >= config$min_mutations
filter_test_regions_by_mut_count <- function(test_region_list, mutations.gr, min_mutations){
  if(min_mutations > 0){
    regions_to_keep <- mutations.gr %>% 
      as_tibble() %>% 
      select(test_region, sampleID) %>% 
      unique() %>% 
      count(test_region) %>% 
      filter(n >= min_mutations) %>%
      pull(test_region)
    
    return(test_region_list[test_region_list %in% regions_to_keep])
  } else {
    return(test_region_list)
  }
}

# Remove overlapping test regions from regions to analyse
filter_overlapping_test_regions <- function(test_region_list, test_regions.gr){
  overlap_diff.df <- test_regions.gr %>%
    findOverlaps(., ., ignore.strand = TRUE) %>% 
    as_tibble() %>% 
    filter(subjectHits != queryHits)
  regions_to_remove <- test_regions.gr[overlap_diff.df$queryHits]$test_region %>% unique()
  
  return(test_region_list[!(test_region_list %in% regions_to_remove)])
}





################################################
# Outcome probability and simulations
################################################

# Function to calculate the likelihood of a cohort outcome given patient gene mutation probabilities
# no_mut_probs - Probabilities of NO mutation for each patient. Slightly faster to supply that from the start than to perform the conversion from mutation probabilites in this function millions of times.
# outcome - Mutated or not (1/0) for each patient
cohort_outcome_prob <- function(no_mut_probs, outcome){
  return(sum(log(abs(outcome - no_mut_probs)))) # sum(log()) instead of log(prod()) (logA + logB = logAB) as the latter can result in too small numbers for R.
}

# Get rank of first value in vector of c(obs_test_value, sim_test_values)
get_obs_rank <- function(obs_test_value, sim_test_values){
  rank(c(obs_test_value, sim_test_values))[1]
}

# Function that gets both no_of_mutated_patients and logp_outcome from a (simulated) cohort outcome
# Used so that we can get get both values from the same set of simulations without saving every outcome vector
sim_recurrence_and_combined_test_values <- function(no_of_patients, pmut){
  sim <- rbinom(no_of_patients, 1, pmut)
  c(sum(sim), cohort_outcome_prob(1 - pmut, sim))
}

# Function to simulate an outcome with a fixed number of mutations
# This uses the standard R sampling method, which is unfortunately really slow when
# replacement=FALSE and the prob weights are used. See improved sim_fixed_mutated_count2 below.
sim_fixed_mutated_count <- function(no_of_patients, no_of_mutations, exp_mut_counts) {
  result <- numeric(no_of_patients)
  result[sample(1:no_of_patients, size=no_of_mutations, prob = exp_mut_counts)] <-  1
  result
}





################################################
# Plotting
################################################

plot_cohort_dist_fit <- function(fit.gamma, obs_val, p_val){
  min_dist <- qgamma(0.01, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"], lower.tail = T)
  max_dist <- qgamma(0.01, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"], lower.tail = F)
  
  if(is.infinite(obs_val)){
    min_lim <- min_dist*0.95
    max_lim <- max_dist*1.05
  } else {
    min_lim <- min(-obs_val, min_dist)*0.95
    max_lim <- max(-obs_val, max_dist)*1.05
  }
  
  plot_max_y <- max(sapply(seq(min_dist, max_dist, by = (max_dist - min_dist) / 100),
                           function(x) dgamma(x, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["rate"])))
  
  p <- tibble(negLogP = c(min_lim, max_lim)) %>% 
    ggplot(aes(x=-obs_val)) +
    stat_function(fun=dgamma, args=list(shape=fit.gamma$estimate["shape"], rate=fit.gamma$estimate["rate"]), aes(fill='simcolor'), geom = 'area', xlim = c(min_lim, max_lim)) +
    xlab('-Log likelihood of cohort mutation distribution')
  if(!is.infinite(obs_val)){
    p <- p + geom_segment(x = -obs_val, y = plot_max_y*0.2, xend = -obs_val, yend = 0, arrow = arrow()) +
      geom_text(x = -obs_val, y = plot_max_y * 0.3, label = paste0('Real cohort\nP = ', formatC(p_val, digits = 2)), hjust="inward")
  } else{
    p <- p + geom_segment(x = max_lim - (max_lim - min_lim)*0.2, y = plot_max_y*0.1, xend = max_lim, yend = plot_max_y*0.1, arrow = arrow()) +
      geom_text(x = max_lim, y = plot_max_y * 0.3, label = paste0('Real cohort\nP = ', formatC(p_val, digits = 2)), hjust="inward")
  }
  p <- p +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(name = '', values = c('simcolor' = '#AEC7E8'), labels = c('simcolor' = 'Simulations')) +
    ggtitle('Cohort mutational skew\ntest result') +
    theme_classic() +
    theme(
      legend.position = c(1, 1),
      legend.justification = c("right", "top"),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y =  element_blank()
      
    )
}


plot_cohort_dist_cumsum <- function(mutprobs_test_region_cohortdist.df, no_of_patients, obs_no_of_mutated_patients, effects_to_keep, no_of_simulations){
  sims <- sapply(1:no_of_simulations, function(x) cumsum(sim_fixed_mutated_count2(no_of_patients, obs_no_of_mutated_patients, mutprobs_test_region_cohortdist.df$exp_mut_count)))
  cumsum_mutated_quantiles <- sapply(1:no_of_patients, function(x) quantile(sims[x,], probs = c(0.05, 0.95))) %>% t() %>% as_tibble() %>% mutate(i = row_number())
  
  label_fun <- function(effect){
    if(effect == 'm'){
      '- Missense'
    } else if(effect == 'n'){
      '- Nonsense'
    } else if(effect == 's'){
      '- Silent'
    } else {
      ''
    }
  }
  
  effect_label <- paste0("Mutated patients (", obs_no_of_mutated_patients, ")")
  effect_label <- ifelse(length(effects_to_keep) == 1 & effects_to_keep[1] == 'na',
                         effect_label,
                         paste0(effect_label, '\n', paste0(c("Included effects:", sapply(effects_to_keep, function(x) label_fun(x))), collapse = '\n')))
  
  obs_cumsum_mutated <- mutprobs_test_region_cohortdist.df %>% arrange(pmut) %>% mutate(i = row_number(), cumsum_mutated = cumsum(mutated))
  
  obs_cumsum_mutated %>%
    mutate(test = 'test') %>% 
    ggplot(aes(fill = test)) +
    geom_ribbon(data = cumsum_mutated_quantiles, aes(x = i, ymin = `5%`, ymax = `95%`, fill="simcolor"), color="#AEC7E8") +
    geom_line(aes(x = i, y = cumsum_mutated, color = "obscolor"))+
    theme_classic() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(name = '', values=c("simcolor"="#AEC7E8"), labels = c("simcolor"="5-95% simulations")) +
    scale_color_manual(name = '', values = c("obscolor" = 'black'), labels = c("obscolor" =  effect_label)) +
    ylab("Mutated patients (cumulative)") +
    xlab('Patients ordered by mutation probability') +
    ggtitle("Cohort mutational skew - cumulative mutated patients")
  # theme(
  #   axis.title.x = element_blank()
  # )
  
}


plot_combined <- function(test_region_arg, cancer_arg, test_regions.gr, mutations_with_effects_for_plotting.df, assembly_arg, annotate_cds_effects, no_of_patients, p.cohort_dist_fit = NULL, p.cohort_dist_cumsum = NULL, highlight_position_min_muts){
  
  no_of_bases <- test_regions.gr %>% filter(test_region == test_region_arg) %>% width() %>% sum()
  if(annotate_cds_effects){
    max_pos <- no_of_bases / 3
    pos_y_lab <- 'Codon'
  } else {
    max_pos <- no_of_bases
    pos_y_lab <- 'Base'
  }
  
  
  plot.df <- mutations_with_effects_for_plotting.df %>% 
    filter(test_region == test_region_arg) %>% 
    left_join(mutprobs_test_region.df %>% arrange(exp_mut_count) %>% mutate(i = row_number()), by = 'sampleID') %>% 
    {if(annotate_cds_effects) mutate(., pos_no = codon_no) else mutate(., pos_no = base_no)}
  
  highlight_positions <- plot.df %>% 
    count(pos_no) %>% 
    filter(n >= highlight_position_min_muts)
  
  if(nrow(highlight_positions) > 0){
    max_pos_for_plot <- 10
    if(nrow(highlight_positions) > max_pos_for_plot){
      highlight_positions_filtered <- highlight_positions %>% 
        arrange(-n) %>% 
        head(n = max_pos_for_plot)
      highlights_filtered_string1 <- paste0(' - top ', max_pos_for_plot, ' darker')
      highlights_filtered_string2 <- paste0(' - top ', max_pos_for_plot)
    } else {
      highlight_positions_filtered <- highlight_positions
      highlights_filtered_string1 <- ''
      highlights_filtered_string2 <- ''
    }
  } else {
    highlights_filtered_string1 <- ''
    highlights_filtered_string2 <- ''
  }
  
  missense_color <- '#f59d47'
  nonsense_color <- '#D62728'
  synonymous_color <- '#4788b5'
  na_color <- 'black'
  
  p <- plot.df %>% 
    ggplot()
  if(nrow(highlight_positions) > 0){
    p <- p +
      geom_hline(data= highlight_positions, aes(yintercept = pos_no), color = 'gray85')
    if(nrow(highlight_positions_filtered) != nrow(highlight_positions)){
      p <- p +
        geom_hline(data= highlight_positions_filtered, aes(yintercept = pos_no), color = 'gray50')
    }
  }
  p <- p +
    geom_point(aes(x=i, y = pos_no, color = effect), size = 1) +
    scale_x_continuous(limits = c(0, no_of_patients), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, max_pos), expand = c(0, 0)) +
    ylab(pos_y_lab) +
    xlab('Patients ordered by mutation probability') +
    ggtitle(paste('Mutation positions', ifelse(nrow(highlight_positions > 0), paste0('(n > ', highlight_position_min_muts - 1,' per ', ifelse(annotate_cds_effects, 'codon', 'base'), ' highlighted',highlights_filtered_string1,')'), ''))) +
    theme_classic() +
    scale_color_manual(name = 'Mutations by effect',
                       values =  c(missense_color, nonsense_color, synonymous_color, na_color),
                       breaks = c('m', 'n', 's', 'na')) +
    theme(legend.position = 'none')
  
  
  if(nrow(highlight_positions) > 0){
    
    if(annotate_cds_effects){
      genome <- load_bsgenome(assembly_arg)
      seq <- get_cds_seq(test_regions.gr, test_region_arg, genome)
      highlighted_pos_info.df <- plot.df %>%
        filter(codon_no %in% highlight_positions_filtered$pos_no) %>%
        count(base_no, codon_no, refnuc, varnuc, strand, effect, pos_no) %>%
        mutate(region_strand = test_regions.gr %>% filter(test_region == test_region_arg) %>% strand() %>% as.character() %>% unique()) %>% 
        mutate(refnuc_coding = ifelse(strand == region_strand, refnuc, str_revcomp(refnuc)),
               varnuc_coding = ifelse(strand == region_strand, varnuc, str_revcomp(varnuc))) %>% 
        rowwise() %>% 
        mutate(label = get_codon_mutation(seq, codon_no, base_no, varnuc_coding)) %>% 
        group_by(label, effect, pos_no) %>%
        summarise(n = sum(n), .groups = 'drop')
    } else {
      highlighted_pos_info.df <- plot.df %>% filter(base_no %in% highlight_positions$pos_no) %>% count(refnuc, base_no, pos_no, varnuc, effect) %>% mutate(label = paste0(refnuc, base_no, varnuc))
    }
    
    p.highlighted_positions <- highlighted_pos_info.df %>% 
      arrange(pos_no, label) %>%
      mutate(label = factor(label, levels = label)) %>% 
      ggplot(aes(x= label, y = n, fill = effect)) +
      geom_col(position = position_dodge2(width = 0.9, preserve = "single")) +
      geom_text(aes(label=n), position=position_dodge(width=0.9), vjust=0.5, hjust = -0.5) +
      theme_classic() +
      ylab("Mutations") +
      xlab(pos_y_lab) +
      ggtitle(paste0('Recurrent mutations\n(n > ', highlight_position_min_muts - 1,' per ', ifelse(annotate_cds_effects, 'codon', 'base'), highlights_filtered_string2, ')')) +
      scale_y_continuous(expand = c(0,0), limits = c(0, max(highlighted_pos_info.df$n)*1.15)) +
      scale_fill_manual(values =  c(missense_color, nonsense_color, synonymous_color, na_color), breaks = c('m', 'n', 's', 'na')) +
      theme(
        legend.position = 'none'
      ) +
      coord_flip()
    

  } else {
    p.highlighted_positions <- NULL
  }
  
  
  effect_count <- plot.df %>% count(effect)
  p.effect_separated_cumsum_muts <- plot.df %>% 
    count(i, effect) %>%
    bind_rows(expand_grid(effect = effect_count$effect,
                          n = 0,
                          i = 0 : no_of_patients)) %>% 
    group_by(i, effect) %>% 
    summarise(n = sum(n), .groups = 'drop') %>% 
    group_by(effect) %>%
    arrange(effect, i) %>%
    mutate(cumsum_muts = cumsum(n), cumsum_mut_ratio = cumsum_muts / max(cumsum_muts)) %>% 
    ggplot(aes(x = i, y = cumsum_mut_ratio, color = effect)) +
    geom_line() +
    theme_classic() +
    ylab('Cumulative portion of mutations') +
    xlab('Patients ordered by mutation probability') +
    ggtitle('Cumulative mutations - separated by effect') +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    scale_color_manual(name = 'Mutations by effect',
                       breaks = c('m', 'n','s', 'na'),
                       labels = c(paste('Mismatch', ifelse('m' %in% effect_count$effect, paste0('(', filter(effect_count, effect == 'm')$n, ')'), '')),
                                  paste('Nonsense', ifelse('n' %in% effect_count$effect, paste0('(', filter(effect_count, effect == 'n')$n, ')'), '')),
                                  paste('Silent', ifelse('s' %in% effect_count$effect, paste0('(', filter(effect_count, effect == 's')$n, ')'), '')),
                                  paste('Mutations', ifelse('na' %in% effect_count$effect, paste0('(', filter(effect_count, effect == 'na')$n, ')'), ''))),
                       values = c(missense_color, nonsense_color, synonymous_color, na_color),
                       guide = guide_legend(override.aes=list(size = 6)))
  
  
  if(!is.null(p.cohort_dist_cumsum) & !is.null(p.cohort_dist_fit)){
    # If cohort dist plots exist, so 2 columns and 3 rows
    p.grid <- plot_grid(
      p.cohort_dist_cumsum,
      p.effect_separated_cumsum_muts,
      p,
      p.cohort_dist_fit,
      NULL,
      p.highlighted_positions,
      nrow = 3,
      align = 'hv',
      axis = 'lrtb',
      byrow = FALSE,
      rel_widths = c(3,2),
      rel_heights = c(2,2,6))
  } else if (nrow(highlight_positions) > 0){
    # If no cohort dist plots, but there is a highlight plot, so still two columns, but only two rows
    p.grid <- plot_grid(
      p.effect_separated_cumsum_muts,
      p,
      NULL,
      p.highlighted_positions,
      nrow = 2,
      align = 'hv',
      axis = 'lrtb',
      byrow = FALSE,
      rel_widths = c(3,2),
      rel_heights = c(2,6))
  } else {
    # If no cohort dist plots, and no highlight plot -> 1 column, 2 rows
    p.grid <- plot_grid(
      p.effect_separated_cumsum_muts,
      p,
      nrow = 2,
      align = 'hv',
      axis = 'lrtb',
      byrow = FALSE,
      rel_heights = c(2,6))
  }
  
  
  title <- ggdraw() + 
    draw_label(
      paste0(test_region_arg ,' (', cancer_arg, ')'),
      fontface = 'bold',
      x = 0.5,
      hjust = 0.5
    )
  plot_grid(
    title, p.grid,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(5, 100)
  )
  
  
  
}





################################################
# Improved sampling method to speed up sampling without replacement, with probability weights.
# https://stackoverflow.com/questions/15113650/faster-weighted-sampling-without-replacement
# Had to put it inside a function that is called in each worker now that I'm using doSNOW. Just exporting faster_sample didn't work (https://stackoverflow.com/a/18245658)
################################################


worker.init <- function(){
  suppressPackageStartupMessages(library(inline))

  src <- 
    '
  int num = as<int>(size), x = as<int>(n);
  Rcpp::NumericVector vx = Rcpp::clone<Rcpp::NumericVector>(x);
  Rcpp::NumericVector pr = Rcpp::clone<Rcpp::NumericVector>(prob);
  Rcpp::NumericVector rnd = rexp(x) / pr;
  for(int i= 0; i<vx.size(); ++i) vx[i] = i;
  std::partial_sort(vx.begin(), vx.begin() + num, vx.end(), Comp(rnd));
  vx = vx[seq(0, num - 1)] + 1;
  return vx;
  '
  incl <- 
    '
  struct Comp{
    Comp(const Rcpp::NumericVector& v ) : _v(v) {}
    bool operator ()(int a, int b) { return _v[a] < _v[b]; }
    const Rcpp::NumericVector& _v;
  };
  '
  assign('faster_sample', cxxfunction(signature(n = "Numeric", size = "integer", prob = "numeric"),
                             src, plugin = "Rcpp", include = incl), .GlobalEnv)

}

sim_fixed_mutated_count2 <- function(no_of_patients, no_of_mutations, exp_mut_counts) {
  result <- numeric(no_of_patients)
  result[faster_sample(no_of_patients, no_of_mutations, exp_mut_counts)] <-  1
  result
}


################################################
# Utility
################################################

str_revcomp <- function(str){
  return(stri_reverse(chartr("TCGA", "AGCT", str)))
}

get_pyr_trinucs <- function(){
  DNA_bases <- c("A", "C", "G", "T")
  pyrs <- c("C", "T")
  pyr_trinucs <- expand.grid(DNA_bases, pyrs, DNA_bases) %>% apply(1, paste0, collapse="")
  return(pyr_trinucs)
}

# Get full range spanned by ranges with the same column_name, e.g. start of first cds to end of last cds for a gene
get_gr_region_range <- function(gr, column_name){
  gr %>%
    as_tibble() %>% 
    group_by((!!sym(column_name)), seqnames, strand) %>% 
    summarise(start = min(start), end = max(end), .groups = "drop") %>% 
    as_granges()
}

# Extend flanks of gr
extend_gr_flanks <- function(gr, flank_upstream, flank_downstream){
  gr %>% 
    {stretch(anchor_5p(.), flank_downstream)} %>% 
    {stretch(anchor_3p(.), flank_upstream)}
}

get_gr_region_range_count <- function(gr, region_column_name, count_column_name){
  gr %>%
    as_tibble() %>%
    group_by((!!sym(region_column_name))) %>%
    summarise((!!sym(count_column_name)) := dplyr::n(), .groups = "drop")
}

# Take a gr, and change seqnames to a combination of seqnames and whatever column "column_name" is.
# Useful for GenomicRanges uperation that automatically reduce when we don't want that.
# For instance, when using setdiff, and reduce is undersired between different genes, this function will create seqnames like "chr1_RPL13A"
# Store the original true chromosome in chrom column
make_fake_paste_column_chrom <- function(gr, column_name, sep = "_"){
  gr %>% 
    as_tibble() %>% 
    mutate(chrom = seqnames) %>%
    mutate(seqnames = paste0(chrom, sep, (!!sym(column_name)))) %>%
    as_granges()
}

# Function to load BSgenome package and return the object
load_bsgenome <- function(assembly){
  if(assembly %in% c('hg19', 'GRCh37')){
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
    return(BSgenome.Hsapiens.UCSC.hg19)
  } else if(assembly %in% c('hg38', 'GRCh38')){
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    return(BSgenome.Hsapiens.UCSC.hg38)
  } else {
    cat('Invalid reference genome. Specify either hg19 or hg38 in config file.\n')
    quit()
  }
}

# Get a GRanges object from BSGenome object genome, optionally only certain chroms.
get_genome_gr <- function(genome, chroms = NA){
  chroms <- ifelse(is.na(chroms), seqnames(chroms), chroms)
  GRanges(seqnames = chroms, strand = '+', ranges = IRanges(start = 1, width = seqlengths(genome)[chroms]))
}

# Get trinuc frequencies on both strands in GRanges object "gr", using BSGenome "genome"
get_genome_trinuc_bgfreqs <- function(genome, gr){
  # Should perhaps stretch out 1 here to get trinucs at mutations right at the end of the regions,
  # but it will make nearly zero difference, and risk a crash if the supplied regions make getSeq try to get something outside of the chromosome boundaries.
  seq <- getSeq(genome, gr) 
  trinuc_bg_counts <- colSums(trinucleotideFrequency(seq) + trinucleotideFrequency(reverseComplement(seq)))
  tibble(trinuc = names(trinuc_bg_counts), count = trinuc_bg_counts)
}

get_chromosome_synonyms <- function(assembly_arg){
  if(assembly_arg %in% c('hg19', 'GRCh37')) assembly_both_formats <- c('hg19', 'GRCh37')
  if(assembly_arg %in% c('hg38', 'GRCh38')) assembly_both_formats <- c('hg38', 'GRCh38')
  return(assembly_both_formats)
}

get_codon_mutation <- function(seq_5_to_3, codon_no, base_no, varnuc_coding_strand){
  codon_table.df <- make_codon_table()
  
  ref_codon <- str_sub(seq_5_to_3, codon_no*3 - 2, codon_no * 3)
  substr(seq_5_to_3, base_no, base_no) <- varnuc_coding_strand
  mut_codon <- str_sub(seq_5_to_3, codon_no*3 - 2, codon_no * 3)
  
  ref_aa <- filter(codon_table.df, codon == ref_codon)$aa
  mut_aa <- filter(codon_table.df, codon == mut_codon)$aa
  
  paste0(ref_aa, codon_no, mut_aa)
}


get_cds_seq <- function(test_regions.gr, test_region_arg, genome){
  region_strand <- test_regions.gr %>% filter(test_region == test_region_arg) %>% strand() %>% as.character() %>% unique()
  test_regions.gr %>%
    filter(test_region == test_region_arg) %>%
    mutate(seq = getSeq(genome, ., as.character = TRUE)) %>%
    as_tibble() %>% 
    { if( region_strand == '+') arrange(., start) else if (region_strand == '-') arrange(., -start) } %>% 
    pull(seq) %>% 
    paste0(collapse = '')
}

# print_logo <- function(){
# logo <- r"{  
#                                                   _--   --_          
#                                                _/           \_       
#                                              /       ___       \     
#                                             /      -     -      \    
#                      ***                   /     /   ***   \     \   
#                     *****                  |    |   *****   |    |   
#                      ***                   \     \   ***   /     /   
#                                             \      -     -      /    
#   ****     *******   ***     ****     **        **   ***      ****   
# ***  ***   *******   ***   ***  ***   ***      ***-__***__--******** 
# ***   **   **        ***   ***   **   ****    ****   ***   ***    ***
#  ****      *****     ***    ****      *****  *****   ***   ***       
#    ****    *****     ***      ****    *** **** ***   ***   ***       
# **   ***   **        ***   **   ***   ***  **  ***   ***   ***    ***
# ***  ***   *******   ***   ***  ***   ***      ***   ***    ******** 
#   ****     *******   ***     ****     ***      ***   ***      ****   
                                                                     
# }"
# cat(logo)
# }

################################################
# Run the main function
################################################

if (!interactive()) {
  main()
}

