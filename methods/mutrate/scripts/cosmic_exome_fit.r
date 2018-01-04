# Fit cosmic signature with exome data using Sigfit R-package

# load libraries
suppressPackageStartupMessages(library(sigfit))
suppressWarnings(library(docopt))
suppressWarnings(library(devtools))

# command line interface (cli)
doc <- "Usage: cosmic_exome_fit.r [-m <mutations>] [-c <cosmic>] [-i <iterations>] [-o <output>]

-m --mutations <mutations> mutational_catalogue.tsv file
-c --cosmic <cosmic> cosmic_signatures.tsv file
-i --iterations <iterations> number of iterations in MCMC
-o --output <output> output folder where sigfit results are kept"

opt <- docopt(doc)

catalogue.file <- opt$mutations
cosmic.exome.file <- opt$cosmic
iterations <- strtoi(opt$iterations)
output.folder <- opt$output


# read mutational catalogue
catalogue.data.frame <- read.csv(file=catalogue.file, colClasses=c("NULL", rep(NA, 96)), 
                                 header=TRUE, sep="\t", row.names = NULL, check.names = FALSE)
catalogue.data.frame[is.na(catalogue.data.frame)] <- 0
mutations.matrix <- data.matrix(catalogue.data.frame)

# read cosmic signatures (in the scope of exome triplet content)
cosmic.exome.data.frame <- read.csv(file=cosmic.exome.file, colClasses=c("NULL", rep(NA, 96)), 
                                    header=TRUE, sep="\t", row.names = NULL, check.names = FALSE)
cosmic.exome.data.frame[is.na(cosmic.exome.data.frame)] <- 0
cosmic.exome.signatures <- data.matrix(cosmic.exome.data.frame)

# MCMC fitting
mcmc_samples_fit <- sigfit::fit_signatures(counts = mutations.matrix, 
                                           signatures = cosmic.exome.signatures, 
                                           iter = iterations, warmup = iterations / 2, 
                                           chains = 1, seed = 1)

# retrieve exposures mean and credible interval from fitting
exposures <- retrieve_pars(mcmc_samples_fit, feature = "exposures", 
                           hpd_prob = 0.90, signature_names = rownames(cosmic.exome.signatures))

# check carefully, as in the first version the method to fill out dataframe cells was wrong!
df.mean <- data.frame(matrix(unlist(exposures$mean), nrow=dim(mutations.matrix)[1], byrow=F))
df.lower <- data.frame(matrix(unlist(exposures$lower_90), nrow=dim(mutations.matrix)[1], byrow=F))
df.upper <- data.frame(matrix(unlist(exposures$upper_90), nrow=dim(mutations.matrix)[1], byrow=F))

# write exposures to .tsv
write.table(df.mean, file=file.path(output.folder, "mean_exposure.tsv"), sep = "\t", quote = FALSE)
write.table(df.lower, file=file.path(output.folder, "lower_exposure.tsv"), sep = "\t", quote = FALSE)
write.table(df.upper, file=file.path(output.folder, "upper_exposure.tsv"), sep = "\t", quote = FALSE)

# plot: summary reconstruction: cohort-wise
sigfit::plot_all(mutations.matrix, 
                 out_path = output.folder, 
                 mcmc_samples = mcmc_samples_fit, 
                 signatures = cosmic.exome.signatures,
                 signature_names = rownames(cosmic.exome.signatures),
                 prefix = "fitting")
