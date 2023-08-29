cat('Installing R dependencies\n')
options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/cran/latest"), timeout=max(1000, getOption("timeout")))

install.packages('BiocManager', version='3.16', dependencies=TRUE)

BiocManager::install("plyranges")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

install.packages('tidyverse', version='1.3.1', dependencies=TRUE)
install.packages('yaml', version='2.3.5')
install.packages('foreach', version='1.5.2')
install.packages('doSNOW', version='1.0.20')
install.packages('stringi', version='1.7.6')
install.packages('inline', version='0.3.19')
install.packages('Rcpp', version='1.0.8.3')
install.packages('fitdistrplus', version='1.1.8')
install.packages('data.table', version='1.14.2')
install.packages('cowplot', version='1.1.1')
install.packages('R.utils', version='2.12.2')
install.packages(c('textshaping','ragg'), dependencies=TRUE)

