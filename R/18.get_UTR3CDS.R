#' Get 3' UTRs and their last CDS regions based on CP sites
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   setup_sqlitedb().
#' @param chr.utr3 An object of [GenomicRanges::GRanges-class], specifying UTR3
#'   GRanges for a chromosome. It must be one element of an output of
#'   [extract_UTR3Anno()]. 
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply. It can be set to NULL or bpparam()
#'
#' @return An object of [GenomicRanges::GRanges-class] containing GRanges for
#'   UTRs with alternative CP sites and the corresponding last CDSs.
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu

get_UTR3CDS <- function(sqlite_db,
                        chr.utr3,
                        BPPARAM = NULL){
    if (missing(sqlite_db) || missing(chr.utr3)){
        stop("sqlite_db and chr.utr3 are required")
    }
    if (length(sqlite_db) != 1 || !file.exists(sqlite_db)){
        stop("The path to the SQLite database is invalid!")
    }
    if (!is(chr.utr3, "GRanges") ||
        !all(chr.utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
        stop("chr.utr3 must be one element of an output of function of",
        "extract_UTR3Anno()")
    }
    
    seqname <- unique(as.character(seqnames(chr.utr3)))
    if (length(seqname) != 1){
        stop("chr.utr3 should contains only one chromosome/scaffold")
    }
    db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_db)
    CPsites <- dbReadTable(db_conn, "CPsites")
    dbDisconnect(db_conn)
    
    chr.CPsites <- CPsites$cpsites_file[CPsites$chr == seqname]
    if (length(chr.CPsites) != 1){
        stop(paste("No CP sites for", seqname))
    }
    
    chr.CPsites <- readRDS(chr.CPsites)
    chr.CPsites <- split(chr.CPsites, names(chr.CPsites))
    if (!is.null(BPPARAM)) {
        utr3.regions <- bptry(bplapply(chr.CPsites,
                                       get_UTR3region,
                                       BPPARAM = BPPARAM))
        while (!all(bpok(utr3.regions))) {
            utr3.regions <- bptry(bplapply(chr.CPsites,
                                           get_UTR3region,
                                           BPREDO = utr3.regions,
                                           BPPARAM = BPPARAM))
        }
    } else {
        utr3.regions <- lapply(chr.CPsites, get_UTR3region)
    }
    utr3.regions <- unlist(GRangesList(utr3.regions))
    
    ## merge utr3.regions with CDS of the same set of transcripts whose 
    ## utr3.regions are considered
    CDS <- chr.utr3 %>% 
        plyranges::filter(feature == "CDS" &
                          transcript %in% unique(utr3.regions$transcript))
    utr3.cds.regions <- c(utr3.regions, CDS)
    utr3.cds.regions
}