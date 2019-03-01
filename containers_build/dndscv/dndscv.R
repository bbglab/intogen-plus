library("dndscv")

args = commandArgs(TRUE)
muts = read.table(gzfile(args[1]), sep = '\t', header = TRUE)

writeOutput = 0
tryCatch({
    result <- dndscv(muts, refdb=file.path(Sys.getenv("INTOGEN_DATASETS"), "dndscv", "RefCDS.rda"))
    writeOutput <- 1
}, error=function(e) {
    message('dndsCV raised an error')
    writeOutput <- 2
})

if (writeOutput == 1) {
    write.table(result$sel_cv, gzfile(args[2]), quote=FALSE, sep='\t', row.names = FALSE)
    write.table(result$annotmuts, gzfile(args[3]), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(result$genemuts, gzfile(args[4]), sep = "\t", quote = FALSE, row.names = FALSE)
} else {
    df = data.frame(
        gene_name=character(),
        n_syn=character(),
        n_mis=character(),
        n_non=character(),
        n_spl=character(),
        n_ind=character(),
        wmis_cv=character(),
        wnon_cv=character(),
        wspl_cv=character(),
        wind_cv=character(),
        pmis_cv=character(),
        ptrunc_cv=character(),
        pallsubs_cv=character(),
        pind_cv=character(),
        qmis_cv=character(),
        qtrunc_cv=character(),
        qallsubs_cv=character(),
        pglobal_cv=character(),
        qglobal_cv=character(),
        stringsAsFactors=FALSE
    )
    write.table(df, gzfile(args[2]), sep = "\t", quote = FALSE, row.names = FALSE)

    df = data.frame(
        sampleID=character(),
        chr=character(),
        pos=character(),
        ref=character(),
        mut=character(),
        gene=character(),
        strand=character(),
        ref_cod=character(),
        mut_cod=character(),
        ref3_cod=character(),
        mut3_cod=character(),
        aachange=character(),
        ntchange=character(),
        codonsub=character(),
        impact=character(),
        pid=character(),
        stringsAsFactors=FALSE
    )
    write.table(df, gzfile(args[3]), sep = "\t", quote = FALSE, row.names = FALSE)

    df = data.frame(
        gene_name=character(),
        n_syn=character(),
        n_mis=character(),
        n_non=character(),
        n_spl=character(),
        exp_syn=character(),
        exp_mis=character(),
        exp_non=character(),
        exp_spl=character(),
        exp_syn_cv=character(),
        stringsAsFactors=FALSE
    )
    write.table(df, gzfile(args[4]), sep = "\t", quote = FALSE, row.names = FALSE)
}

