library("dndscv")

args = commandArgs(TRUE)
muts = read.table(gzfile(args[1]), sep = '\t', header = TRUE)
target_genes = read.table(args[2], sep = '\t', header = FALSE)

# TODO capture only the specific errors we know about
writeOutput = 0
tryCatch({
    # TODO add the RefCDS.rda directly to the image instead of using the datasest ?
    # or pass it explicitly
    #NORMAL execution
    #result <- dndscv(muts, refdb=file.path(Sys.getenv("INTOGEN_DATASETS"), "dndscv", "RefCDS.rda"))
    #PANEL execution
    result <- dndscv(muts, gene_list=target_genes, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf)
    # if including refdb as a parameter, the gene_list is not found in the DB. Check format before using it. It is required for panels?
    #result <- dndscv(muts, gene_list=target_genes, max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, refdb=file.path(Sys.getenv("INTOGEN_DATASETS"), "dndscv", "RefCDS.rda"))
    writeOutput <- 1
}, error=function(e) {
    message('dndsCV raised an error')
    writeOutput <- 2
})

if (writeOutput == 1) {
    write.table(result$sel_cv, gzfile(args[3]), quote=FALSE, sep='\t', row.names = FALSE)
    write.table(result$annotmuts, gzfile(args[4]), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(result$genemuts, gzfile(args[5]), sep = "\t", quote = FALSE, row.names = FALSE)
} else {
    # TODO have the empty files existing and just copy them
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
    write.table(df, gzfile(args[3]), sep = "\t", quote = FALSE, row.names = FALSE)

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
    write.table(df, gzfile(args[4]), sep = "\t", quote = FALSE, row.names = FALSE)

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
    write.table(df, gzfile(args[5]), sep = "\t", quote = FALSE, row.names = FALSE)
}

