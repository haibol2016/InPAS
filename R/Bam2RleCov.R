#' Get read coverage in the Rle format from a BAM file
#'
#' @param bamfile A character(1) file name of the 'BAM' file to be processed.
#' @param index A character(1) name of the index file of the 'BAM' file being
#'   processed; this is given without the '.bai' extension
#' @param tag A character(1) name tags used to label the input bamfile
#' @param param An instance of [Rsamtools::ScanBamParam-class]. This influences 
#'   what fields and which records are imported. Use of which requires that a 
#'   BAM index file (filename.bai) exists.
#' @param mapqFilter An integer(1), mapping quality filtering threshold. 
#'   Default is 255, good for uniquely mapping reads by RNA aligner, STAR.
#' @param run.type A character(1), sequencing read layout: single end or 
#'   paired end reads?
#' @param removeScaffolds A logical(1), whether the scaffolds should be removed
#'   from the genome
#' @param hugeData A logical(1). Is this dataset consume too much memory? if it is
#'   TRUE, the coverage will be saved into tempfiles.
#' @param outdir A directory with write permission for storing intermediate
#'   coverage data. This is required for huge data, i.e., when hugeData is set
#'   to TRUE
#'
#' @return If hugeData is set "FALSE", the function returns a list containing 
#'   read coverage and total depth for a bedgraph file. The first list element
#'   contains a list of Rle instances of [S4Vectors::Rle-class] representing
#'   read coverage for each chromosome, with the chromosome names starting 
#'   with "chr" as names. The second element contains total depth. Otherwise,
#'   the returned value is a character(1), with tag as the name.
#'   
#' @export
#' 
#' @importFrom Rsamtools scanBamFlag ScanBamParam
#' @importClassesFrom Rsamtools ScanBamParam
#' @importFrom GenomicAlignments readGAlignments readGAlignmentPairs coverage
#'

Bam2RleCov <- function(bamfile, 
                       index = bamfile,
                       tag,
                       param = NULL,
                       mapqFilter = 255,
                       run.type = 
                           c("single_end", "paired_end"),
                       removeScaffolds = FALSE,
                       hugeData = TRUE,
                       outdir = NULL){
    if (missing(bamfile) || length(bamfile) != 1){
        stop("A single bamfile is required")
    }

    if (missing(tag) || length(tag) != 1){
        stop("A single name tag for the bamfile is required")
    }
    
    if (file.exists(bamfile)){
        stop("The bamfile does NOT exist!")
    }
    
    if (hugeData && is.null(outdir))
    {
        stop("A writable output directory is", 
             " required for store huge coverage data!")
    }
    
    run.type <- match.arg(run.type)
    
    bf <- BamFile(bamfile)
    
    ## get chromosome names and their length from the BAM header
    bam_header <- scanBamHeader(bf, what = "targets")
    bam_header <-  bam_header$targets
    
    if (removeScaffolds) {
        bam_header <- bam_header[grepl("^(chr)?(\\d+|[XY])$", 
                                     names(bam_header))]
    }
    
    flag <- scanBamFlag(isPaired = NA, 
                        isProperPair = NA, 
                        isUnmappedQuery = FALSE,  
                        hasUnmappedMate = NA, 
                        isMinusStrand = NA, 
                        isMateMinusStrand = NA, 
                        isFirstMateRead = NA, 
                        isSecondMateRead = NA, 
                        isSecondaryAlignment = FALSE,
                        isNotPassingQualityControls = FALSE,
                        isDuplicate = NA)
    # RleList <- list(tag = list())
    # 
    # if (missing(param)){
    #     for (chr in names(bam_header)) {
    #         param <- ScanBamParam(flag = flag, 
    #                               simpleCigar = FALSE, 
    #                               reverseComplement = FALSE, 
    #                               tag = character(0), 
    #                               tagFilter = list(), 
    #                               what = character(0), 
    #                               which = 
    #                                   GRanges(seqnames = chr,
    #                                     ranges = IRanges(1, bam_header[chr])),
    #                               mapqFilter = mapqFilter) 
    #         
    #         if (run.type == "single_end")
    #         {
    #             gal <- readGAlignments(bf, param=param)
    #         } else {
    #             gal <- readGAlignmentPairs(bf, param=param)
    #         }
    #         seqlevels(gal) <- chr
    #         cov <- coverage(gal)
    #         RleList$tag[[chr]] <- cov
    #     }
    # } else {
    #     for (chr in names(bam_header)) {
    #         bamWhich(param) <- GRanges(seqnames = chr,
    #                                    ranges = IRanges(1, bam_header[chr]))
    #         if (run.type == "single_end")
    #         {
    #             gal <- readGAlignments(bf, param=param)
    #         } else {
    #             gal <- readGAlignmentPairs(bf, param=param)
    #         }
    #         seqlevels(gal) <- chr
    #         cov <- unlist(coverage(gal))
    #         RleList$tag[[chr]] <- cov
    #     }
    # }
    # names(RleList) <- tag

    RleList <- lapply(names(bam_header), function(.x){
        if (is.null(param)) {
            param <- ScanBamParam(flag = flag, 
                                  simpleCigar = FALSE, 
                                  reverseComplement = FALSE, 
                                  tag = character(0), 
                                  tagFilter = list(), 
                                  what = character(0), 
                                  which = 
                                      GRanges(seqnames = .x,
                                 ranges = IRanges(1, bam_header[.x])),
                                  mapqFilter = mapqFilter) 
            if (run.type == "single_end")
            {
                gal <- readGAlignments(bf, param=param)
            } else {
                gal <- readGAlignmentPairs(bf, param=param)
            }
        } else {
            bamWhich(param) <- GRanges(seqnames = chr,
                                    ranges = IRanges(1, bam_header[chr]))
            if (run.type == "single_end"){
                gal <- readGAlignments(bf, param=param)
            } else {
                gal <- readGAlignmentPairs(bf, param=param)
            }
        }
        seqlevels(gal) <- .x
        cov <- coverage(gal)
        cov
    })
    cov_list <- lapply(RleList, function(.x){
       .x[[1]]
    })
    names(cov_list) <- names(bam_header)
    
    ## total depth per sample
    depth <- sapply(cov_list, function(.cvg) {
        sum(as.double(runValue(.cvg)) * runLength(.cvg))
    })
    depth <- sum(depth)
    names(depth) <- tag
    
    RleList <- list(cov = cov_list, depth = depth)
    names(RleList)[1] <- tag
    
    if (hugeData)
    {
        ## save coverage as files into a directory named by chromosome ID
        filenames <- lapply(names(bam_header), function(chr){
            od <- file.path(outdir, chr)
            if (!dir.exists(od))
            {
                dir.create(od, recursive = TRUE)
            }
            od <- normalizePath(od)
            filename <- file.path(od, 
                                  paste(tag, chr, 
                                        Sys.Date(), 
                                        "RleCov.RDS",sep = "_"))
            saveRDS(RleList[[tag]][chr], file = filename)
            filename
        })
        names(filenames) <- names(bam_header)
        filenames <- list(fn = filenames, depth = depth)
        names(filenames)[1] <- tag
        
        # return files in a list of list, first element is chromosome 
        # ID-named list containing filenames, second element contains 
        # total coverage depth of the sample
        filenames
    } else {
        RleList
    }
}





