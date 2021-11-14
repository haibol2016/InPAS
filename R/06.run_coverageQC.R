
#' Quality control on read coverage over gene bodies and 3UTRs
#'
#' Calculate coverage over gene bodies and 3UTRs. This function is used for
#' quality control of the coverage.The coverage rate can show the complexity of
#' RNA-seq library.
#'
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param TxDb An object of [GenomicFeatures::TxDb-class]
#' @param edb An object of [ensembldb::EnsDb-class]
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param removeScaffolds A logical(1) vector, whether the scaffolds should be
#'   removed from the genome If you use a TxDb containing alternative
#'   scaffolds, you'd better to remove the scaffolds.
#' @param cutoff_readsNum cutoff reads number. If the coverage in the location
#'   is greater than cutoff_readsNum,the location will be treated as covered by
#'   signal
#' @param cutoff_expdGene_cvgRate cutoff_expdGene_cvgRate and
#'   cutoff_expdGene_sampleRate are  the  parameters used to calculate which
#'   gene is expressed in all input dataset. cutoff_expdGene_cvgRateset the
#'   cutoff value for the coverage rate of each gene;
#'   cutoff_expdGene_sampleRateset the cutoff value for ratio of numbers of
#'   expressed and all samples for each gene. for example, by default,
#'   cutoff_expdGene_cvgRate=0.1 and cutoff_expdGene_sampleRate=0.5,suppose
#'   there are 4 samples, for one gene, if the coverage rates by base are:0.05,
#'   0.12, 0.2, 0.17, this gene will be count as expressed gene because
#'   mean(c(0.05,0.12, 0.2, 0.17) > cutoff_expdGene_cvgRate) >
#'   cutoff_expdGene_sampleRate if the coverage rates by base are: 0.05, 0.12,
#'   0.07, 0.17, this gene will be count as un-expressed gene because
#'   mean(c(0.05, 0.12, 0.07, 0.17) > cutoff_expdGene_cvgRate)<=
#'   cutoff_expdGene_sampleRate
#' @param cutoff_expdGene_sampleRate See cutoff_expdGene_cvgRate
#' @param which an object of [GenomicRanges::GRanges-class] or NULL. If it is
#'   not NULL, only the exons overlapping the given ranges are used. For fast 
#'   data quality control, set which to Granges for one or a few large 
#'   chromosomes.
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply. It can be set to NULL or bpparam() 
#' @param ... Not used yet
#'
#' @return A data frame with colnames: gene.coverage.rate: coverage per base for
#'   all genes, expressed.gene.coverage.rate: coverage per base for expressed
#'   genes, UTR3.coverage.rate: coverage per base for all 3'
#'   UTRs,UTR3.expressed.gene.subset.coverage. rate: coverage per base for 3'
#'   UTRs of expressed genes. and rownames: the names of coverage
#' @import GenomicRanges
#' @importFrom BiocParallel bptry bpok bplapply bpparam
#'
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' if (interactive()) {
#'    library("BSgenome.Mmusculus.UCSC.mm10")
#'    library("TxDb.Mmusculus.UCSC.mm10.knownGene")
#'    library("EnsDb.Mmusculus.v79")
#'    
#'    genome <- BSgenome.Mmusculus.UCSC.mm10
#'    TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#'    edb <- EnsDb.Mmusculus.v79
#'    
#'    bedgraphs <- system.file("extdata",c("Baf3.extract.bedgraph",
#'                                         "UM15.extract.bedgraph"), 
#'                            package = "InPAS")
#'    tags <- c("Baf3", "UM15")
#'    metadata <- data.frame(tag = tags, 
#'                           condition = c("Baf3", "UM15"),
#'                           bedgraph_file = bedgraphs)
#'    outdir = tempdir()
#'    write.table(metadata, file =file.path(outdir, "metadata.txt"), 
#'                sep = "\t", quote = FALSE, row.names = FALSE)
#'    
#'    sqlite_db <- setup_sqlitedb(metadata = file.path(outdir, 
#'                                                     "metadata.txt"),
#'                                outdir)
#'    coverage <- list()
#'    for (i in seq_along(bedgraphs)){
#'    coverage[[tags[i]]] <- get_ssRleCov(bedgraph = bedgraphs[i],
#'                             tag = tags[i],
#'                             genome = genome,
#'                             sqlite_db = sqlite_db,
#'                             outdir = outdir,
#'                             removeScaffolds = TRUE,
#'                             BPPARAM = NULL)
#'    }
#'    coverage_files <- assemble_allCov(sqlite_db, 
#'                                     outdir, 
#'                                     genome, 
#'                                     removeScaffolds = FALSE)
#'    run_coverageQC(sqlite_db, TxDb, edb, genome,
#'                   removeScaffolds = TRUE,
#'                   which = GRanges("chr6",
#'                      ranges = IRanges(98013000, 140678000)))
#' }

run_coverageQC <- function(sqlite_db, TxDb, edb, genome,
                         cutoff_readsNum = 1,
                         cutoff_expdGene_cvgRate = 0.1,
                         cutoff_expdGene_sampleRate = 0.5,
                         removeScaffolds = FALSE,
                         BPPARAM = NULL,
                         which = NULL, ...) {
  stopifnot(is(TxDb, "TxDb"))
  stopifnot(is(edb, "EnsDb"))
  stopifnot(is(genome, "BSgenome"))
  if (missing(sqlite_db) || length(sqlite_db) != 1 ||
      !file.exists(sqlite_db)) {
    stop("sqlite_db, a path to the SQLite database is required!")
  }

  seqnames <- trim_seqnames(genome, removeScaffolds)
  seqLen <- get_seqLen(genome, removeScaffolds)
  
  if (length(which) > 0) {
    stopifnot(is(which, "GRanges"))
    message("strand information will be ignored.")
    query_seqname <- unique(as.character(seqnames(which)))
    seqnames <- seqnames[seqnames %in% query_seqname]
    if(length(seqnames) < 1){
      stop("The seqnames defined by 'which' is not valid")
    }
  }

  db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
  coverage_files <- dbGetQuery(db_conn, 
                               "SELECT * FROM sample_coverage")
  coverage_files <- coverage_files[coverage_files$chr %in% seqnames, ]
  dbDisconnect(db_conn)
  
  if (nrow(coverage_files) < 1){
    stop("No coverage files in the SQLite database!")
  }
  
  coverage <- tags <- unique(coverage_files$tag)
  names(coverage) <- tags
  
  seqRle <- sapply(seqLen, function(.ele) {Rle(0, .ele)},
                   simplify = FALSE)

  tx <- parse_TxDb(TxDb, edb, removeScaffolds = removeScaffolds)
  exon <- tx %>% plyranges::select(exon, feature, gene)
  UTR3 <- tx %>%
    plyranges::filter(feature %in% c("lastutr3", "utr3")) %>%
    plyranges::mutate(feature = "utr3")

  if (length(which) > 0) {
    ol <- GenomicRanges::findOverlaps(exon, which, ignore.strand = TRUE)
    tx <- exon[sort(unique(queryHits(ol)))]
    ol <- GenomicRanges::findOverlaps(UTR3, which, ignore.strand = TRUE)
    UTR3 <- UTR3[sort(unique(queryHits(ol)))]
  }

  ## reduce exon and UTR3 for each gene
  reduce_by_gene <- function(gr) {
    gr <- gr %>%
      plyranges::group_by(seqnames, gene) %>%
      plyranges::reduce_ranges_directed()
    gr
  }

  exon <- reduce_by_gene(exon)
  UTR3 <- reduce_by_gene(UTR3)

  feature.cvg <- function(cov, feature, seqnames, seqRle) {
    ## read in sample coverage data
    sample_cov <- list()
    for (chr in seqnames){
      f <- coverage_files$tag == cov & coverage_files$chr == chr
      if (sum(f) == 1)
      {
        sample_cov[[chr]] <- 
          readRDS(coverage_files$coverage_file[f])
      } else { ## no coverage
        sample_cov[[chr]] <- seqRle[[chr]]
      }
    }
    cov <- sample_cov
    cov <- cov[seqnames]
    names(cov) <- seqnames
    cvg.base <- mapply(function(cov.seq, feature.seq) {
      vw <- Views(cov.seq,
        start = start(feature.seq),
        end = end(feature.seq)
      )
      viewApply(vw, function(.ele) {
        sum(unlist(.ele > cutoff_readsNum))
      })
    }, cov[seqnames], feature[seqnames], SIMPLIFY = FALSE)
    stopifnot(identical(
      sapply(cvg.base, length),
      sapply(feature, length)))
    unlist(cvg.base, use.names = FALSE)
  }

  exon.s <- split(exon, as.character(seqnames(exon)))
  seqnames <- seqnames[seqnames %in% names(exon.s)]
  if (length(seqnames) < 1) {
    stop("No chromosome is selected.")
  }
  exon.s <- exon.s[seqnames]
  seqRle <- seqRle[seqnames]
  exon.unlist <-
    if (!is(exon.s, "GRangesList")) GRangesList(exon.s) else exon.s
  exon.unlist <- unlist(exon.s)

  if (is.null(BPPARAM)){
    exon.cvg <- lapply(coverage, feature.cvg,
                       feature = exon.s,
                       seqnames = seqnames, seqRle = seqRle)
  } else {
    exon.cvg <-bptry(bplapply(coverage, feature.cvg,
                       feature = exon.s,
                       seqnames = seqnames, seqRle = seqRle,
                       BPPARAM = BPPARAM))
    while (!all(bpok(exon.cvg))) {
      exon.cvg <-bptry(bplapply(coverage, feature.cvg,
                                feature = exon.s,
                                seqnames = seqnames, seqRle = seqRle, 
                                BPREDO = exon.cvg, BPPARAM = BPPARAM))
    }
  }


  names(exon.cvg) <- names(coverage)

  ## coverage per sample in each column
  exon.cvg <- do.call(cbind, exon.cvg)
  exon.cvg.wid <- cbind(exon.cvg, gene.width = width(exon.unlist))

  ## summary of exonic coverage and width per genes
  gene.cvg <- rowsum(exon.cvg.wid, exon.unlist$gene, reorder = FALSE)
  gene.expd.rate <- gene.cvg[, -ncol(gene.cvg)] / gene.cvg[, ncol(gene.cvg)]
  gene.expd.rate <- gene.expd.rate > cutoff_expdGene_cvgRate

  gene.expd.overall <- rowSums(gene.expd.rate) >
    (cutoff_expdGene_sampleRate * ncol(gene.expd.rate))
  gene.expd.overall.rate <- mean(gene.expd.overall)

  ## overall averaged gene.cvg.rate per sample
  gene.cvg.rate <- colSums(gene.cvg[, -ncol(gene.cvg)]) / 
    sum(gene.cvg[, ncol(gene.cvg)])
  if (gene.expd.overall.rate > 0) {
    gene.cvg.rate.subsetBy.expd.gene <-
      colSums(gene.cvg[gene.expd.overall, -ncol(gene.cvg)]) /
        sum(gene.cvg[gene.expd.overall, ncol(gene.cvg)])
  } else {
    gene.cvg.rate.subsetBy.expd.gene <- rep(0, length(coverage))
    names(gene.cvg.rate.subsetBy.expd.gene) <- names(coverage)
  }

  UTR3.s <- split(UTR3, as.character(seqnames(UTR3)))
  seqnames <- seqnames[seqnames %in% names(UTR3.s)]
  if (length(seqnames) < 1) {
    stop("No chromosome is selected.")
  }
  UTR3.s <- UTR3.s[seqnames]
  seqRle <- seqRle[seqnames]
  UTR3.unlist <- if (!is(UTR3.s, "GRangesList")) GRangesList(UTR3.s) else UTR3.s
  UTR3.unlist <- unlist(UTR3.s)

  UTR3.cvg <- lapply(coverage, feature.cvg,
    feature = UTR3.s,
    seqnames = seqnames, seqRle = seqRle
  )
  names(UTR3.cvg) <- names(coverage)
  UTR3.cvg <- do.call(cbind, UTR3.cvg)
  UTR3.cvg.rate <- colSums(UTR3.cvg) / sum(width(UTR3.unlist))
  idx <- UTR3.unlist$gene %in% names(gene.expd.overall[gene.expd.overall])
  if (sum(idx) > 0) {
    UTR3.cvg.rate.subsetBy.expd.gene <-
      colSums(UTR3.cvg[idx, , drop = FALSE]) / sum(width(UTR3.unlist[idx]))
  } else {
    UTR3.cvg.rate.subsetBy.expd.gene <- rep(0, length(coverage))
    names(UTR3.cvg.rate.subsetBy.expd.gene) <- names(coverage)
  }

  data.frame(gene.coverage.rate = gene.cvg.rate,
             expressed.gene.coverage.rate = 
                          gene.cvg.rate.subsetBy.expd.gene,
             UTR3.coverage.rate = UTR3.cvg.rate,
             UTR3.expressed.gene.subset.coverage.rate =
                          UTR3.cvg.rate.subsetBy.expd.gene)
}
