#' extract coverage from a bedgraph file
#'
#' extract coverage from a bedgraph file
#'
#' @param bedgraph A bedGraph file name with its accessible path
#' @param genome an object [BSgenome::BSgenome-class]
#' @param removeScaffolds logical(1), whether the scaffolds should be removed
#'   from the genome
#' @param hugeData is this dataset consume too much memory? if it is TRUE, the
#'   coverage will be saved into tempfiles.
#' @param outdir a directory with write permission for storing intermediate
#'   coverage data. This is required for huge data, i.e., when hugeData is set
#'   to TRUE
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply. It can be set to NULL or bpparam()
#'   
#' @return an object of [S4Vectors::Rle-class] for a sample coverage
#' @importFrom dplyr as_tibble mutate filter arrange bind_rows group_by
#'   left_join summarise n bind_cols syms desc pull
#' @importFrom GenomeInfoDb seqlengths seqlevelsStyle mapSeqlevels
#' @import readr
#' @import BiocParallel
#' @import S4Vectors GenomicRanges
#' @importFrom BSgenome getSeq matchPWM
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @examples 
#'   library(BSgenome.Mmusculus.UCSC.mm10)
#'   bedgraph <- system.file("extdata",
#'     "Baf3.extract.bedgraph",
#'     package = "InPAS"
#'   )
#'   genome <- BSgenome.Mmusculus.UCSC.mm10
#'   coverage <- getSingleSampleCov(bedgraph,
#'     genome,
#'     removeScaffolds = TRUE,
#'     hugeData = FALSE,
#'     outdir = tempdir()
#'   )
#'   coverage <- getSingleSampleCov(bedgraph,
#'     genome,
#'     removeScaffolds = TRUE,
#'     hugeData = TRUE,
#'     BPPARAM = BiocParallel::bpparam(),
#'     outdir = tempdir()
#'   )

getSingleSampleCov <- function(bedgraph, 
                   genome, 
                   removeScaffolds = FALSE,
                   hugeData = TRUE,
                   BPPARAM = NULL,
                   outdir = NULL) {
    if (missing(genome)) {
        stop("genome is required.")
    }
    if (!is(genome, "BSgenome")) {
        stop("genome must be an object of BSgenome.")
    }
    if (length(bedgraph) != 1) {
        stop("length of bedgraph must be 1")
    }
    if (!file.exists(bedgraph)) {
        stop("bedgraph file does not exist")
    }
    if (hugeData && is.null(outdir))
    {
        stop("A writable output directory is", 
             " required for store huge coverage data!")
    }

    seqnames.bedfile <-
        as.data.frame(read_tsv(
            bedgraph,
            comment = "#",
            col_names = FALSE, skip = 0,
            col_types = cols(
                X1 = col_factor(),
                X2 = "-",
                X3 = "-",
                X4 = "-"
            )
        ))[, 1]
    
    seqnames <- trimSeqnames(genome, removeScaffolds)
    
    seqStyle <- seqlevelsStyle(genome)
    seqStyle.bed <- seqlevelsStyle(levels(seqnames.bedfile))
    
    if (any(seqStyle.bed != seqStyle)) {
        stop("seqlevelsStyle of genome is different from bedgraph file.")
     }
    seqnames <- sort(intersect(levels(seqnames.bedfile), seqnames))
    if (length(seqnames) < 1) {
        stop(paste(
            "there is no intersect seqname in",
            bedgraph, "and genome"
        ))
    }
    
    seqLen <- seqLen(genome, removeScaffolds)
    
    summaryFunction <- function(seqname) {
        seqL <- seqLen[seqname]
        
        ## apply some trick here by using "FALSE"
        lines2read <- Rle(c(FALSE, seqnames.bedfile == seqname))
        true <- which(runValue(lines2read))
        skip <- runLength(lines2read)[true - 1]
        skip[1] <- skip[1] - 1
        nrow <- runLength(lines2read)[true]
        
        dat <- read_tsv(bedgraph,
                        comment = "#",
                        col_names = FALSE,
                        skip = skip[1],
                        n_max = nrow[1],
                        col_types = cols(
                            X1 = "-",
                            X2 = "i",
                            X3 = "i",
                            X4 = "d"
                        )
        )
        
        if (length(true) > 1) {
            cul_skip <- skip[1] + nrow[1]
            for (i in 2:length(true)) {
                cul_skip <- cul_skip + skip[i]
                lines <-
                    read_tsv(bedgraph,
                             comment = "#",
                             col_names = FALSE,
                             skip = cul_skip,
                             n_max = nrow[i],
                             col_types = cols(
                                 X1 = "-",
                                 X2 = "i",
                                 X3 = "i",
                                 X4 = "d"
                             )
                    )
                dat <- dplyr::bind_rows(dat, lines)
                cul_skip <- cul_skip + nrow[i]
            }
        }
        
        # convert bedgraph 0-based index to 1-based index for GRanges
        dat <- dat %>% dplyr::mutate(X2 = X2 + 1)
        
        ## convert uncovered regions as gaps
        ## This step can be avoided when generate bedgraph with bedtools:
        ## bedtools genomecov -bga -split -ibam  <coordinate-sorted BAM file>
        gaps <-
            as_tibble(as.data.frame(gaps(IRanges(
                dplyr::pull(dat, 1),
                dplyr::pull(dat, 2)
            ),
            start = 1, end = seqL
            )))
        if (nrow(gaps) > 0) {
            gaps <- dplyr::bind_cols(gaps[, 1:2],
                                     V4 = rep(0, nrow(gaps))
            )
            colnames(gaps) <- colnames(dat)
            dat <- dplyr::bind_rows(dat, gaps)
        }
        dat <- dat %>%
            dplyr::arrange(X2) %>%
            dplyr::mutate(width = X3 - X2 + 1)
        Rle(dplyr::pull(dat, 3), dplyr::pull(dat, width))
    }
    
    if (!is.null(BPPARAM))
    {
        cvg <- bptry(bplapply(seqnames, summaryFunction, BPPARAM = BPPARAM))
        while (!all(bpok(cvg))) {
            cvg <- bptry(bplapply(seqnames, summaryFunction,
            BPREDO = cvg, BPPARAM = BPPARAM))
        }
    } else {
        cvg <- lapply(seqnames, summaryFunction)
    }
    names(cvg) <- seqnames
    if (hugeData)
    {
        ## save files into a directory called "coverage"
        if (!dir.exists(outdir))
        {
            dir.create(outdir, recursive = TRUE)
        }
        
        filename <- file.path(normalizePath(outdir), 
                             paste(basename(bedgraph), basename(tempfile()), 
                                    "coverage.RData",sep = "_"))
        save(list = "cvg", file = filename)
        filename 
    } else {
        cvg
    }
}

#' Combine coverage files
#' 
#' Combine per-sample coverage files of an experiment into a list
#'
#' @param coverageFiles paths to coverage files, absolute paths or path relative
#'   to the current working directory 
#' @param tags unique sample name tags for each corresponding bedGraph file
#'
#' @return A list of coverage with tags as names
#' @export


getCovFileList <- function(coverageFiles, tags){
    if (missing(tags) || missing(coverageFiles)) {
        stop("tags and bedgraphs are required.")
    }
    if (length(tags) != length(coverageFiles)) {
        stop("lengths of tags and bedgraphs must be identical.")
    }
    if (!is.character(tags)) {
        stop("tags must be a character vector")
    }
    if (any(duplicated(tags))) {
        stop("tags must be unique")
    }
    if (any(!file.exists(coverageFiles))) {
        stop("Not all bedgraphs exist")
    }
    
    covFileList <- as.list(coverageFiles)
    names(covFileList) <- tags
    covFileList
}


#' get coverage for small-scale experiments
#'
#' @param bedgraphs a vector of character, containing paths to bedGraph files,
#'   absolute paths or path relative to the current working directory 
#' @param tags unique sample name tags matching each bedGraph file
#' @param genome an object [BSgenome::BSgenome-class]
#' @param removeScaffolds logical(1), whether the scaffolds should be removed
#'   from the genome
#' @param hugeData a vector of logical(1), indicating whether the data is huge.
#'   For experiments with more than 3 samples, hugeData = TRUE should be 
#'   considered.
#' @param outdir a directory with write permission for storing intermediate
#'   coverage data. This is required for hughe data, i.e., when hugeData is set
#'   to TRUE
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply. It can be set to NULL or bpparam()
#'
#' @return a list of coverage files or objects of [S4Vectors::Rle-class] 
#' @export

getCov4SmallExperiment <- function(bedgraphs, 
                             tags,
                             genome, 
                             removeScaffolds = FALSE,
                             hugeData = FALSE,
                             outdir = NULL,
                             BPPARAM = NULL){
    if (missing(tags) || missing(bedgraphs)) {
        stop("tags and bedgraphs are required.")
    }
    if (length(tags) != length(bedgraphs)) {
        stop("length of tags and bedgraphs should be identical.")
    }
    if (!is.character(tags)) {
        stop("tags must be a character vector")
    }
    if (any(duplicated(tags))) {
        stop("tags must be unique")
    }
    if (any(!file.exists(bedgraphs))) {
        stop("Not all bedgraphs exist")
    }
    if (missing(genome)) {
        stop("genome is required.")
    }
    if (!is(genome, "BSgenome")) {
        stop("genome must be an object of BSgenome.")
    }
    if (length(bedgraphs) >= 4){
        stop("For an experiment with more than 3 samples, please use the ",
        "functions getSingleSampleCov() and getCovFileList subsequentially")
    }
    if (hugeData && is.null(outdir))
    {
        stop("A writable output directory is", 
             " required for store huge coverage data!")
    }
    
    cov <- lapply(bedgraphs, function(.x){
        getSingleSampleCov(bedgraph = .x, 
                           genome = genome, 
                           removeScaffolds = removeScaffolds,
                           hugeData = hugeData,
                           BPPARAM = BPPARAM,
                           outdir = outdir)
    })
    names(cov) <- tags
    cov
}