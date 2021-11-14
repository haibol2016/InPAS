#' Estimate the CP sites for UTRs on a given chromosome
#'
#' Estimate the CP sites for UTRs on a given chromosome
#' 
#' @param seqname A character(1) vector, specifying a chromososome/scaffold name
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param utr3 An object of [GenomicRanges::GRanges-class]. Output of [extract_UTR3Anno()] for a chromosome/scaffold
#' @param background A character(1) vector, the range for calculating cutoff
#'   threshold of local background. It can be "same_as_long_coverage_threshold",
#'   "1K", "5K","10K", or "50K".
#' @param z2s one element of an output of [setup_CPsSearch()] for Z-score cutoff
#'   values, which is the output of [get_zScoreCutoff()]
#' @param depth.weight A named vector. One element of an output of
#'   [setup_CPsSearch()] for coverage depth weight, which is the output of
#'   [get_depthWeight()]
#' @param genome A [BSgenome::BSgenome-class] object
#' @param MINSIZE A integer(1) vector, specifying the minimal length in bp of a
#'   short/proximal 3' UTR. Default, 10
#' @param window_size An integer(1) vector, the window size for novel distal or 
#'   proximal CP site searching. default: 100.
#' @param search_point_START A integer(1) vector, starting point relative to 
#'   the 5' extremity of 3' UTRs for searching for proximal CP sites 
#' @param search_point_END A integer(1) vector, ending point relative to the 3'
#'   extremity of 3' UTRs for searching for proximal CP sites 
#' @param cutStart An integer(1) vector a numeric(1) vector. What percentage or
#'   how many nucleotides should be removed from 5' extremities before searching
#'   for CP sites? It can be a decimal between 0, and 1, or an integer greater 
#'   than 1. 0.1 means 10 percent, 25 means cut first 25 bases
#' @param cutEnd An integer(1) vector a numeric(1) vector. What percentage or
#'   how many nucleotides should be removed from 5' extremities before searching
#'   for CP sites? It can be a decimal between 0, and 1, or an integer greater 
#'   than 1. 0.1 means 10 percent, 25 means cut first 25 bases
#' @param adjust_distal_polyA_end A logical(1) vector. If true, distal CP sites
#'   are subject to adjustment by the Naive Bayes classifier from the
#'   [cleanUpdTSeq::cleanUpdTSeq-package]
#' @param coverage_threshold An integer(1) vector, specifying the cutoff 
#'   threshold of coverage for first 100 nucleotides. If the coverage of first 
#'   100 nucleotides is lower than coverage_threshold, that transcript will be
#'   not considered for further analysis. Default, 5.
#' @param long_coverage_threshold An integer(1) vector, specifying the cutoff 
#'   threshold of coverage for the terminal of long form 3' UTRs. If the coverage
#'   of first 100 nucleotides is lower than coverage_threshold, that transcript
#'   will be not considered for further analysis. Default, 2.
#' @param PolyA_PWM  An R object for a position weight matrix (PWM) for a hexamer
#'   polyadenylation signal (PAS), such as AAUAAA. 
#' @param classifier An R object for Naive Bayes classifier model, like the one
#'   in the cleanUpdTSeq package.
#' @param classifier_cutoff A numeric(1) vector. A cutoff of probability that a
#'   site is classified as true CP sites. The value should be between 0.5 and 1.
#'   Default, 0.8.
#' @param shift_range An integer(1) vector, specifying a shift range for 
#'   adjusting the proximal and distal CP sites. Default, 100. It determines the
#'   range flanking the candidate CP sites to search the most likely real 
#'   CP sites.
#' @param step An integer (1) vector, specifying the step size used for adjusting
#'   the proximal or distal CP sites using the Naive Bayes classifier from the
#'   cleanUpdTSeq package. Default 1. It can be in the range of 1 to 10.
#' @param two_way A logical (1), indicating whether the proximal CP sites are
#'   searched from both directions or not.
#' @param hugeData A logical(1), indicating whether it is huge data
#' @param outdir A character(1) vector, a path with write permission for storing 
#'   the CP sites. If it doesn't exist, it will be created.
#' @param silence logical(1), indicating whether progress is reported or not. By
#'   default, FALSE
#'
#' @return An object of [GenomicRanges::GRanges-class] containing distal and
#'   proximal CP site information for each 3' UTR
#' @seealso [search_proximalCPs()], [adjust_proximalCPs()],
#'   [adjust_proximalCPsByPWM()], [adjust_proximalCPsByNBC()],
#'   [get_PAscore()], [get_PAscore2()]
#' @import GenomicRanges
#' @importFrom BSgenome getSeq matchPWM
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples 
#' if (interactive()) {
#'    library(BSgenome.Mmusculus.UCSC.mm10)
#'    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#'    genome <- BSgenome.Mmusculus.UCSC.mm10
#'    TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#'    
#'    ## load UTR3 annotation and convert it into a GRangesList
#'    data(utr3.mm10)
#'    utr3 <- split(utr3.mm10, seqnames(utr3.mm10))
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
#'                                          "metadata.txt"), outdir)
#'    coverage <- list()
#'    for (i in seq_along(bedgraphs)) {
#'    coverage[[tags[i]]] <- get_ssRleCov(bedgraph = bedgraphs[i],
#'                             tag = tags[i],
#'                             genome = genome,
#'                             sqlite_db = sqlite_db,
#'                             outdir = outdir,
#'                             removeScaffolds = TRUE,
#'                             BPPARAM = NULL)}
#'    coverage_files <- assemble_allCov(sqlite_db, 
#'                                     outdir, 
#'                                     genome, 
#'                                     removeScaffolds = TRUE)
#'    data4CPsSearch <- setup_CPsSearch(sqlite_db,
#'                                      genome,
#'                                      utr3,
#'                                      background = "10K",
#'                                      TxDb = TxDb,
#'                                      removeScaffolds = TRUE,
#'                                      BPPARAM = NULL,
#'                                      hugeData = TRUE,
#'                                      outdir = outdir)
#'    ## polyA_PWM
#'    load(system.file("extdata", "polyA.rda", package = "InPAS"))
#'    
#'    ## load the Naive Bayes classifier model from the cleanUpdTSeq package
#'    library(cleanUpdTSeq)
#'    data(classifier)
#'    
#'    CPs <- search_CPs(seqname = "chr6",
#'                      sqlite_db = sqlite_db, 
#'                      utr3 = utr3,
#'                      background = data4CPsSearch$background, 
#'                      z2s = data4CPsSearch$z2s,
#'                      depth.weight = data4CPsSearch$depth.weight,
#'                      genome = genome, 
#'                      MINSIZE = 10, 
#'                      window_size = 100,
#'                      search_point_START =50,
#'                      search_point_END = NA,
#'                      cutStart = 10, 
#'                      cutEnd = 0,
#'                      adjust_distal_polyA_end = TRUE,
#'                      coverage_threshold = 5,
#'                      long_coverage_threshold = 2,
#'                      PolyA_PWM = pwm, 
#'                      classifier = classifier,
#'                      classifier_cutoff = 0.8,
#'                      shift_range = 100,
#'                      step = 5,
#'                      two_way = FALSE,
#'                      hugeData = TRUE,
#'                      outdir = outdir)
#' }
#' 

search_CPs <- function(seqname,
                       sqlite_db, 
                       utr3,
                       background, 
                       z2s,
                       depth.weight,
                       genome, 
                       MINSIZE = 10, 
                       window_size = 100,
                       search_point_START = 50,
                       search_point_END = NA,
                       cutStart = 10, 
                       cutEnd = 0,
                       adjust_distal_polyA_end = TRUE,
                       coverage_threshold = 5,
                       long_coverage_threshold = 2,
                       PolyA_PWM = NA, 
                       classifier = NA,
                       classifier_cutoff = 0.8,
                       shift_range = window_size,
                       step = 1,
                       two_way = FALSE,
                       hugeData = TRUE,
                       outdir, 
                       silence = FALSE) {

  if (!is.na(PolyA_PWM)[1]) {
    if (!is(PolyA_PWM, "matrix")) stop("PolyA_PWM must be matrix")
    if (any(rownames(PolyA_PWM) != c("A", "C", "G", "T"))) {
      stop("rownames of PolyA_PWM must be c('A', 'C', 'G', 'T')")
    }
  }
  if (missing(genome) || missing(utr3)) {
    stop("genome and utr3 are required.")
  }
  if (!is(genome, "BSgenome")) {
    stop("genome must be an object of BSgenome.")
  }
  if (!is(utr3, "GRangesList") |
      !all(utr3[[1]]$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
    stop("utr3 must be output of function of extract_UTR3Anno()")
  }
  
  if (seqlevelsStyle(utr3[[1]]) != seqlevelsStyle(genome)) {
    stop("the seqlevelsStyle of utr3 must be same as genome")
  }
  
  if (missing(outdir)){
    stop("An explicit output directory is required")
  } else {
    outdir <- file.path(outdir, "CPsites.out")
    if (!dir.exists(outdir)){
      dir.create(outdir, recursive = TRUE, 
                 showWarnings = FALSE)
    }
    outdir <- normalizePath(outdir)
  }
  
  db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
  utr3_coverage <- dbReadTable(db_conn, "utr3_total_coverage")
  dbDisconnect(db_conn)
  chr.cov <- utr3_coverage$coverage_file[utr3_coverage$chr == seqname]
  
  if (length(chr.cov)!= 1){
    stop("The seqname is not included in the UTR3TotalCov")
  }
  
  ## load coverage file for chr
  chr.cov <- readRDS(file = chr.cov)
  if(!is.list(chr.cov)){
    stop("Something wrong when loading big data. ", 
         "Maybe the tempfile is broken!")
  }

  chr.cov <- chr.cov[sapply(chr.cov, mean) > 0]
  if (length(chr.cov) == 0) {
    return(NULL)
  }
  
  ## utr3 is GRangesList; subset by seqname to get utr3 GRanges with feature ==
  ## "CDS"
  curr_UTR.gr <- utr3 <- utr3[[seqname]] %>%
    plyranges::filter(feature != "CDS")
  if (length(utr3) == 0){return(NULL)}

  utr3.utr <- curr_UTR.gr[curr_UTR.gr$feature == "utr3"]
  utr3.gap <- curr_UTR.gr[curr_UTR.gr$feature == "next.exon.gap"]
  co <- countOverlaps(utr3.gap, utr3.utr, maxgap = 1, 
                      ignore.strand = TRUE)
  utr3.gap <- utr3.gap[co > 1]
  curr_UTR.gr$conn_next_utr3 <-
    curr_UTR.gr$transcript %in% utr3.gap$transcript
  
  ## remove utr3 Granges without coverage
  curr_UTR.gr <- curr_UTR.gr[names(curr_UTR.gr) %in% names(chr.cov)]
  if (length(curr_UTR.gr) == 0){
    return(NULL)
  }
  
  curr_UTR <- split(curr_UTR.gr, curr_UTR.gr$transcript)
  conn_next_utr3 <- sapply(curr_UTR, function(.UTR) {
    .UTR$conn_next_utr3[1]
  })
  
  chr.cov.merge <- lapply(curr_UTR, function(.UTR) {
    .UTR <- .UTR[order(start(.UTR))]
      chr.utr3TotalCov <- chr.cov[names(.UTR)]
      chr.utr3TotalCov <-
        mapply(function(.covList, .start, .end, .property) {
          # set names for each position
          .posList <- .start:.end
          
          # if not a matrix
          if (length(dim(.covList)) == 0 && !is.null(.covList)) {
            .covList <- t(.covList)
          } 
          rownames(.covList) <- paste(.property, .posList, sep = "_SEP_")
          .covList
        }, chr.utr3TotalCov, start(.UTR), end(.UTR), 
        .UTR$feature, SIMPLIFY = FALSE)
      
      chr.utr3TotalCov <- do.call(rbind, chr.utr3TotalCov)
      ## reverse the negative strand
      if (as.character(strand(.UTR))[1] == "-") { 
        chr.utr3TotalCov <-
          chr.utr3TotalCov[rev(rownames(chr.utr3TotalCov)), , drop = FALSE]
      }
      
      # if the range of "cutstart" is (0, 1), percentage; 
      # otherwise, absolute bases
      if (!is.na(cutStart)) {
        if (cutStart < 1) {
          cutStart <- floor(length(chr.utr3TotalCov) * cutStart)
        }
        if (cutStart > 0) {chr.utr3TotalCov <-
                                chr.utr3TotalCov[-(1:cutStart), , drop = FALSE]
        }
      }
      chr.utr3TotalCov
  })
  ## chr.cov.merge should be a list with named numeric
  if (!silence) message("chromsome ", seqname, " coverage merged.\n")

  ## filter UTR3 based on coverage of the first 100 bases
  coverage_quality <- sapply(chr.cov.merge, function(.ele) {
    if (nrow(.ele[grepl("utr3_SEP_", rownames(.ele)), ,
      drop = FALSE]) > MINSIZE) {
      any(colMeans(.ele[1:min(nrow(.ele), 100), ,
        drop = FALSE]) > coverage_threshold)
    } else {
      FALSE
    }})

  chr.cov.merge <- chr.cov.merge[coverage_quality]
  conn_next_utr3 <- conn_next_utr3[coverage_quality]
  if (!silence) {
    message("chromsome ", seqname, 
            " quality filtered by first 100nt coverage.\n")
  }

  if (length(chr.cov.merge) > 0) {
    ## Step 1: search distal CP sites
    if (!silence) message("chromsome ", seqname, " distal search ... start.\n")
    curr_UTR <- curr_UTR[names(chr.cov.merge)]
    chr.abun <- search_distalCPs(chr.cov.merge, conn_next_utr3, curr_UTR,
                                 window_size, depth.weight,
                                 long_coverage_threshold,
                                 background, z2s)
    if (!silence) {message("chromsome ", seqname, " distal search ... done.\n")}

    ## Step 2: adjust distal CP sites
    if (!silence) message("chromsome ", seqname, " distal adjust ... start.\n")
    if (adjust_distal_polyA_end && is(classifier, "PASclassifier")) {
      chr.abun <- adjust_distalCPs(chr.abun, classifier,
                                   classifier_cutoff,
                                   shift_range, genome, step)
      if (!silence) { message("chromsome ", seqname, 
                              " distal adjust ... done.\n")}
    }

    ## Step 3: search proximal CP sites
    if (!silence) message("chromsome ", seqname, 
                          " proximal search ... start.\n")
    chr.abun <- search_proximalCPs(chr.abun, curr_UTR,
                                   window_size, MINSIZE,
                                   cutEnd, search_point_START,
                                   search_point_END, two_way)
    if (!silence) {message("chromsome ", seqname, 
                           " proximal searched ... done.\n")}

    ## Step 4: adjust proximal CP sites
    if (!silence) message("chromsome ", seqname, 
                          " proximal adjust ... start.\n")
    if (is(PolyA_PWM, "matrix") || is(classifier, "PASclassifier")) {
      chr.abun <- adjust_proximalCPs(chr.abun, MINSIZE,
                                     PolyA_PWM, genome,
                                     classifier, classifier_cutoff,
                                     shift_range, search_point_START,
                                     step)
      if (!silence) {message("chromsome ", seqname, 
                             " proximal adjust ... done.\n")}
    }
    chr.abun <- polish_CPs(chr.abun)
  } else {
    chr.abun <- NULL
  }
  
  if (!is.null(chr.abun)){
    utr3.shorten.UTR <- utr3 %>%
      plyranges::filter(feature == "utr3") %>%
      plyranges::select(-c(feature)) %>%
      plyranges::filter(transcript %in% rownames(chr.abun))
    chr.abun <- chr.abun[utr3.shorten.UTR$transcript, , drop = FALSE]
    
    utr3.shorten.UTR$fit_value <- unlist(chr.abun[, "fit_value"])
    utr3.shorten.UTR$Predicted_Proximal_APA <-
      unlist(chr.abun[, "Predicted_Proximal_APA"])
    utr3.shorten.UTR$Predicted_Distal_APA <-
      unlist(chr.abun[, "Predicted_Distal_APA"])
    utr3.shorten.UTR$type <- unlist(chr.abun[, "type"])
    utr3.shorten.UTR$utr3start <- unlist(chr.abun[, "utr3start"])
    utr3.shorten.UTR$utr3end <- unlist(chr.abun[, "utr3end"])
    utr3.shorten.UTR <-
      utr3.shorten.UTR[!is.na(utr3.shorten.UTR$Predicted_Proximal_APA)]
  } else {
    utr3.shorten.UTR <- NULL
  }
  
  ## save chromosome-wise CP sites
  if (!is.null(utr3.shorten.UTR)) {
    filename <- file.path(outdir, paste0(seqname, "_CPsites.RDS"))
    saveRDS(utr3.shorten.UTR, file = filename)
    
    insist_execute_sqlite <- fix_dbLockError()
    tryCatch({
      db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
      
      # remove existing record fpr this seqname
      insist_execute_sqlite(db_conn, "BEGIN IMMEDIATE TRANSACTION")
      # Rollback on failure
      on.exit(try(dbExecute(db_conn, "ROLLBACK TRANSACTION")))
      res <- dbSendStatement(db_conn, 
                             paste0("DELETE FROM   
                                CPsites WHERE chr = '", seqname, "';"))
      dbClearResult(res)
      dbExecute(db_conn, "COMMIT TRANSACTION")
      
      # insert new data
      insist_execute_sqlite(db_conn, "BEGIN IMMEDIATE TRANSACTION")
      # Rollback on failure
      on.exit(try(dbExecute(db_conn, "ROLLBACK TRANSACTION")))
      res <- dbSendStatement(db_conn, 
                             paste0("INSERT INTO 
                                CPsites (chr, cpsites_file) 
                                VALUES ('", seqname, "','", filename,"');"))
      dbClearResult(res)
      dbExecute(db_conn, "COMMIT TRANSACTION")
      
      # Don't Roll back on success
      on.exit(NULL)
    }, error = function(e) {
      print(paste(conditionMessage(e)))
    }, finally = dbDisconnect(db_conn))
  }
  utr3.shorten.UTR
}