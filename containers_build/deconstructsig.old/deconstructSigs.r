#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

## load libraries
library(deconstructSigs)
library(gplots)

# read our input file
ttype=args[1]
x<-read.table(paste(ttype, sep=""), header=T, sep="\t", colClasses=c(NA,NA,NA,NA,NA,"NULL"))
colnames(x) <- c("chr", "pos", "ref", "alt", "Sample")
x$chr <- as.factor(x$chr)
head(x)

# Select genome reference
bsg = NULL
if (Sys.getenv("INTOGEN_GENOME") == "hg38") {
    bsg = BSgenome.Hsapiens.UCSC.hg38::Hsapiens
}

if (Sys.getenv("INTOGEN_GENOME") == "hg19") {
    bsg = BSgenome.Hsapiens.UCSC.hg19::Hsapiens
}

# Convert to deconstructSigs input
# this step generates the matrix suitable for the program
sigs.input <- mut.to.sigs.input(mut.ref = x, 
                                sample.id = "Sample", 
                                chr = "chr", 
                                pos = "pos", 
                                ref = "ref", 
                                alt = "alt",
                                bsg = bsg)
# now we run deconstructSigs for each sample in our input list
flag = 0
for (sample in unique(x$Sample))
{

    if (nrow(x[which(x$Sample==sample),]) > 25) # here I restrict the analysis for only samples with more than 25 mutations
    {
        test = whichSignatures(tumor.ref = sigs.input, 
                               signatures.ref = signatures.cosmic, 
                               sample.id = sample,
                               contexts.needed = TRUE,
                               tri.counts.method = 'exome2genome')

        a = test$weights # save the weights for each signature. 
        a['SSE']  = round(sqrt(sum(test$diff * test$diff)), digits = 3) # compute the error rate
        a['mutation_count'] = nrow(x[which(x$Sample==sample),]) # number of mutations
        # append the results of each sample in to dataframe
        if (flag == 0){total = a; flag=1}
        else{total <- rbind(total, a)}
    }
}
if (flag==1){
# prepare heatmap
new = as.matrix(total[ , grepl( "Signature" , names( total ) ) ])
pdf(paste("heatmap_deconstruct", ".pdf", sep=""))
heatmap.2(new,dendrogram='none', Rowv=TRUE, Colv=TRUE,trace='none', cexRow=0.5, cexCol=0.5)
dev.off()

# prepare CSV file
myDF <- cbind(sample_id = rownames(total), total) # assign row names
rownames(myDF) <- NULL
# write the output to a file
write.table(myDF, file=paste("signatures_weight", ".csv", sep=""), sep="\t", col.names = TRUE, row.names=FALSE)
}
