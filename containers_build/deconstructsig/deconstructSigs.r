#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# load libraries

library(deconstructSigs)
suppressMessages(library(gplots))

# read our input file

ttype=args[1]
x<-read.table(paste(ttype, sep=""), header=T, sep="\t", colClasses=c(NA,NA,NA,NA,NA,"NULL"))
colnames(x) <- c("chr", "pos", "ref", "alt", "Sample")
x$chr <- as.factor(x$chr)

# read genome build as 2nd argument

build=args[2]

# read output file as 3rd argument

output=args[3]

# select genome reference

bsg = NULL

if (build == "hg38") 
{
    bsg = BSgenome.Hsapiens.UCSC.hg38::Hsapiens
}

if (build == "hg19") 
{
    bsg = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
}

message("Creating mutational catalogue...")

# convert to deconstructSigs input
# this step generates the matrix suitable for the program

suppressWarnings(
sigs.input <- mut.to.sigs.input(mut.ref = x, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = bsg)
)

# if a sample has less than completionThreshold mutations, 
# then add the difference spread by the expected number of 
# mutations per channel of the global profile

globalProfile <- colSums(sigs.input)
globalProfile <- globalProfile / sum(globalProfile)
completionThreshold <- 100

for (sample in row.names(sigs.input))
{
    w <- sigs.input[sample,]
    n <- sum(w)
    if (n < completionThreshold)
    {
        k <- completionThreshold - n
        sigs.input[sample,] <- sigs.input[sample,] + (k * globalProfile)
    }
}

message("Running signature fitting...")

# run deconstructSigs for each sample in our input list

flag = 0
for (sample in unique(x$Sample))
{
    if (nrow(x[which(x$Sample==sample),]) > 0)
    {
        test = whichSignatures(tumor.ref = sigs.input, 
                               signatures.ref = signatures.cosmic, 
                               sample.id = sample,
                               signature.cutoff = 0.05,
                               contexts.needed = TRUE,
                               tri.counts.method = 'exome2genome')

        a = test$weights # save the weights for each signature. 
        a['SSE']  = round(sqrt(sum(test$diff * test$diff)), digits = 3) # compute the error rate
        a['mutation_count'] = nrow(x[which(x$Sample==sample),]) # number of mutations
        
        # append the results of each sample to dataframe

        if (flag == 0){total = a; flag=1}
        else{total <- rbind(total, a)}
    }
}

if (flag==1){

message("Plotting heatmap...")

# prepare heatmap

new = as.matrix(total[ , grepl( "Signature" , names( total ) ) ])
pdf(paste("heatmap_deconstruct", ".pdf", sep=""))
heatmap.2(new, dendrogram='none', Rowv=TRUE, Colv=TRUE,trace='none', cexRow=0.5, cexCol=0.5)
dev.off()

message("Saving dataframe...")

# prepare CSV file

myDF <- cbind(sample_id = rownames(total), total)
rownames(myDF) <- NULL

# write the output to a file

write.table(myDF, file=output, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

}
