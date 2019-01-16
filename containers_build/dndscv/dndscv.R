library("dndscv")

args = commandArgs(TRUE)
muts = read.table(gzfile(args[1]), sep = '\t', header = TRUE)
result = dndscv(muts, refdb=file.path(Sys.getenv("INTOGEN_DATASETS"), "dndscv", "RefCDS.rda"))
write.table(result$sel_cv, gzfile(args[2]), quote=FALSE, sep='\t', row.names = FALSE)
write.table(result$annotmuts, gzfile(args[3]), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(result$genemuts, gzfile(args[4]), sep = "\t", quote = FALSE, row.names = FALSE)

