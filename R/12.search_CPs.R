#' Estimate the CP sites for UTRs on a given chromosome
#'
#' Estimate the CP sites for UTRs on a given chromosome
#'
#' @param seqname A character(1) vector, specifying a chromososome/scaffold name
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   [setup_sqlitedb()].
#' @param chr.utr3 An object of [GenomicRanges::GRanges-class]. An element of
#'   the output of [extract_UTR3Anno()].
#' @param genome A [BSgenome::BSgenome-class] object
#' @param MINSIZE A integer(1) vector, specifying the minimal length in bp of a
#'   short/proximal 3' UTR. Default, 10
#' @param window_size An integer(1) vector, the window size for novel distal or
#'   proximal CP site searching. default: 200.
#' @param search_point_START A integer(1) vector, starting point relative to
#'   the 5' extremity of 3' UTRs for searching for proximal CP sites
#' @param search_point_END A integer(1) vector, ending point relative to the 3'
#'   extremity of 3' UTRs for searching for proximal CP sites
#' @param cutEnd An integer(1) vector a numeric(1) vector. What percentage or
#'   how many nucleotides should be removed from 5' extremities before searching
#'   for proximal CP sites? It can be a decimal between 0, and 1, or an integer
#'   greater than 1. 0.1 means 10 percent, 25 means cut first 25 bases
#' @param filter.last A logical(1), whether to filter out the last valley, which
#'   is likely the 3' end of the longer 3' UTR if no novel distal CP site is
#'   detected and the 3' end excluded by setting cutEnd/search_point_END is small.
#' @param adjust_distal_polyA_end A logical(1) vector. If true, distal CP sites
#'   are subject to adjustment by the Naive Bayes classifier from the
#'   [cleanUpdTSeq::cleanUpdTSeq-package]
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
#'   adjusting the proximal and distal CP sites. Default, 50. It determines the
#'   range flanking the candidate CP sites to search the most likely real
#'   CP sites.
#' @param step An integer (1) vector, specifying the step size used for adjusting
#'   the proximal or distal CP sites using the Naive Bayes classifier from the
#'   cleanUpdTSeq package. Default 1. It can be in the range of 1 to 10.
#' @param outdir A character(1) vector, a path with write permission for storing
#'   the CP sites. If it doesn't exist, it will be created.
#' @param future.chunk.size The average number of elements per future
#'   ("chunk"). If Inf, then all elements are processed in a single future.
#'   If NULL, then argument future.scheduling = 1 is used by default. Users can
#'   set future.chunk.size = total number of elements/number of cores set for
#'   the backend. See the future.apply package for details. Default, 50. This
#'   parameter is used to split the candidate 3' UTRs for alternative SP sites
#'   search.
#' @param silence A logical(1), indicating whether progress is reported or not.
#'   By default, FALSE
#' @param mc.cores An integer(1), number of cores for making multicore clusters
#'   or socket clusters using \pkg{batchtools}, and for [parallel::mclapply()]
#' @param cluster_type A character (1) vector, indicating the type of cluster
#'   job management systems. Options are "interactive","multicore", "torque",
#'   "slurm", "sge", "lsf", "openlava", and "socket". see
#'   \href{https://mllg.github.io/batchtools/articles/batchtools.html}{batchtools vignette}
#' @param template_file A charcter(1) vector, indicating the template file for
#'   job submitting scripts when cluster_type is set to "torque", "slurm",
#'   "sge", "lsf", or "openlava".
#' @param resources A named list specifying the computing resources when
#'   cluster_type is set to "torque", "slurm", "sge", "lsf", or "openlava". See
#'   \href{https://mllg.github.io/batchtools/articles/batchtools.html}{batchtools vignette}
#' @return An object of [GenomicRanges::GRanges-class] containing distal and
#'   proximal CP site information for each 3' UTR
#' @seealso [search_proximalCPs()], [adjust_proximalCPs()],
#'   [adjust_proximalCPsByPWM()], [adjust_proximalCPsByNBC()],
#'   [get_PAscore()], [get_PAscore2()]
#' @import GenomicRanges
#' @importFrom batchtools makeRegistry batchMap submitJobs waitForJobs
#'   reduceResultsList makeClusterFunctionsInteractive makeClusterFunctionsLSF
#'   makeClusterFunctionsMulticore makeClusterFunctionsOpenLava
#'   makeClusterFunctionsSGE makeClusterFunctionsSlurm
#'   makeClusterFunctionsSocket makeClusterFunctionsTORQUE
#' @importFrom BSgenome getSeq matchPWM
#' @importFrom parallelly availableCores
#' @importFrom future.apply future_mapply future_lapply
#' @export
#' @author Jianhong Ou, Haibo Liu
#'
#' @examples
#' if (interactive()) {
#'   library(BSgenome.Mmusculus.UCSC.mm10)
#'   library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#'   genome <- BSgenome.Mmusculus.UCSC.mm10
#'   TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
#'
#'   ## load UTR3 annotation and convert it into a GRangesList
#'   data(utr3.mm10)
#'   utr3 <- split(utr3.mm10, seqnames(utr3.mm10), drop = TRUE)
#'
#'   bedgraphs <- system.file("extdata", c(
#'     "Baf3.extract.bedgraph",
#'     "UM15.extract.bedgraph"
#'   ),
#'   package = "InPAS"
#'   )
#'   tags <- c("Baf3", "UM15")
#'   metadata <- data.frame(
#'     tag = tags,
#'     condition = c("Baf3", "UM15"),
#'     bedgraph_file = bedgraphs
#'   )
#'   outdir <- tempdir()
#'   write.table(metadata,
#'     file = file.path(outdir, "metadata.txt"),
#'     sep = "\t", quote = FALSE, row.names = FALSE
#'   )
#'
#'   sqlite_db <- setup_sqlitedb(metadata = file.path(
#'     outdir,
#'     "metadata.txt"
#'   ), outdir)
#'   addLockName(filename = tempfile())
#'   coverage <- list()
#'   for (i in seq_along(bedgraphs)) {
#'     coverage[[tags[i]]] <- get_ssRleCov(
#'       bedgraph = bedgraphs[i],
#'       tag = tags[i],
#'       genome = genome,
#'       sqlite_db = sqlite_db,
#'       outdir = outdir,
#'       chr2exclude = "chrM"
#'     )
#'   }
#'   data4CPsSearch <- setup_CPsSearch(sqlite_db,
#'     genome,
#'     chr.utr3 = utr3[["chr6"]],
#'     seqname = "chr6",
#'     background = "10K",
#'     TxDb = TxDb,
#'     hugeData = TRUE,
#'     outdir = outdir,
#'     minZ = 2,
#'     cutStart = 10,
#'     MINSIZE = 10,
#'     coverage_threshold = 5
#'   )
#'   ## polyA_PWM
#'   load(system.file("extdata", "polyA.rda", package = "InPAS"))
#'
#'   ## load the Naive Bayes classifier model from the cleanUpdTSeq package
#'   library(cleanUpdTSeq)
#'   data(classifier)
#'   ## the following setting just for demo.
#'   if (.Platform$OS.type == "window") {
#'     plan(multisession)
#'   } else {
#'     plan(multicore)
#'   }
#'   CPs <- search_CPs(
#'     seqname = "chr6",
#'     sqlite_db = sqlite_db,
#'     chr.utr3 = utr3[["chr6"]],
#'     genome = genome,
#'     MINSIZE = 10,
#'     window_size = 100,
#'     search_point_START = 50,
#'     search_point_END = NA,
#'     cutEnd = 0,
#'     filter.last = TRUE,
#'     adjust_distal_polyA_end = TRUE,
#'     long_coverage_threshold = 2,
#'     PolyA_PWM = pwm,
#'     classifier = classifier,
#'     classifier_cutoff = 0.8,
#'     shift_range = 100,
#'     step = 5,
#'     outdir = outdir
#'   )
#' }
search_CPs <- function(seqname,
                       sqlite_db,
                       chr.utr3,
                       genome = getInPASGenome(),
                       MINSIZE = 10,
                       window_size = 200,
                       search_point_START = 100,
                       search_point_END = NA,
                       cutEnd = 200,
                       filter.last = TRUE,
                       adjust_distal_polyA_end = TRUE,
                       long_coverage_threshold = 2,
                       PolyA_PWM = NA,
                       classifier = NA,
                       classifier_cutoff = 0.8,
                       shift_range = 50,
                       step = 1,
                       outdir = getInPASOutputDirectory(),
                       silence = FALSE,
                       cluster_type = c(
                         "interactive", "multicore",
                         "torque", "slurm", "sge",
                         "lsf", "openlava", "socket"
                       ),
                       template_file = NULL,
                       mc.cores = 1,
                       future.chunk.size = 50,
                       resources = list(
                         walltime = 3600 * 8, num.cpu = 4,
                         mpp = 1024 * 4, queue = "long"
                       )) {
  if (!is.null(mc.cores) && !is.numeric(mc.cores)) {
    stop("mc.cores must be an integer(1)")
  } else if (.Platform$OS.type == "windows") {
    mc.cores <- 1
  } else if (is.null(mc.cores)) {
    mc.cores <- availableCores() - 1
  } else {
    mc.cores <- min(availableCores() - 1, mc.cores)
  }

  if (!is.na(PolyA_PWM)[1]) {
    if (!is(PolyA_PWM, "matrix")) stop("PolyA_PWM must be matrix")
    if (any(rownames(PolyA_PWM) != c("A", "C", "G", "T"))) {
      stop("rownames of PolyA_PWM must be c('A', 'C', 'G', 'T')")
    }
  }

  if (!is.na(PolyA_PWM) && !is.na(classifier)) {
    stop("PolyA_PWM and classifier can't be set at the same time for adjusting!")
  }

  if (missing(genome) || missing(chr.utr3)) {
    stop("genome and chr.utr3 are required.")
  }
  if (!is(genome, "BSgenome")) {
    stop("genome must be an object of BSgenome.")
  }
  if (!is(chr.utr3, "GRanges") |
    !all(chr.utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
    stop("utr3 must be output of function of extract_UTR3Anno()")
  }

  if (seqlevelsStyle(chr.utr3) != seqlevelsStyle(genome)) {
    stop("the seqlevelsStyle of utr3 must be same as genome")
  }

  if (missing(outdir)) {
    stop("An explicit output directory is required")
  } else {
    outdir <- file.path(outdir, "006.CPsites.out")
    if (!dir.exists(outdir)) {
      dir.create(outdir,
        recursive = TRUE,
        showWarnings = TRUE
      )
    }
    outdir <- normalizePath(outdir)
  }
  lock_filename <- getLockName()
  if (!file.exists(lock_filename)) {
    stop(
      "lock_filename must be an existing file.",
      "Please call addLockName() first!"
    )
  }
  file_lock <- flock::lock(lock_filename)
  db_conn <- dbConnect(
    drv = RSQLite::SQLite(),
    dbname = sqlite_db
  )
  utr3_coverage <- dbReadTable(db_conn, "utr3_total_coverage")
  dbDisconnect(db_conn)
  flock::unlock(file_lock)

  chr.cov <- utr3_coverage$coverage_file[utr3_coverage$chr == seqname]

  if (length(chr.cov) != 1) {
    stop("The seqname is not included in the UTR3TotalCov")
  }

  ## load CPsSearch_data file, "seqname_data4CPsSearch.RDS", for chr = seqname
  chr.cov <- readRDS(file = chr.cov)
  if (!is.list(chr.cov)) {
    stop(
      "Something wrong when loading big data.\n",
      "Maybe the data prepared for CP site search is broken!"
    )
  }

  curr_UTR <- chr.cov$chr.utr3
  chr.cov.merge <- chr.cov$chr.cov.merge
  conn_next_utr3 <- chr.cov$conn_next_utr3
  depth.weight <- chr.cov$depth.weight
  background <- chr.cov$background
  z2s <- chr.cov$z2s

  get_CPsites <- function(seqname,
                          chr.cov.merge,
                          conn_next_utr3,
                          curr_UTR,
                          depth.weight,
                          background,
                          z2s,
                          window_size,
                          long_coverage_threshold,
                          adjust_distal_polyA_end,
                          classifier,
                          classifier_cutoff,
                          shift_range,
                          genome,
                          step,
                          MINSIZE,
                          search_point_START,
                          search_point_END,
                          cutEnd,
                          filter.last,
                          PolyA_PWM,
                          outdir,
                          silence) {
    ## Step 1: search distal CP sites
    if (!silence) {
      message(
        "chromsome ", seqname,
        " distal search start at ", date(), ".\n"
      )
    }
    chr.abun <- search_distalCPs(chr.cov.merge,
      conn_next_utr3,
      curr_UTR,
      window_size = window_size,
      depth.weight,
      long_coverage_threshold =
        long_coverage_threshold,
      background,
      z2s
    )

    if (!silence) {
      message(
        "chromsome ", seqname, " distal search done at ",
        date(), ".\n"
      )
    }

    ## Step 2: adjust distal CP sites
    if (!silence) {
      message(
        "chromsome ", seqname, " distal adjust start at ",
        date(), ".\n"
      )
    }
    if (adjust_distal_polyA_end) {
      chr.abun <- adjust_distalCPs(
        chr.abun,
        classifier,
        classifier_cutoff,
        shift_range,
        genome, step
      )
    }
    if (!silence) {
      message(
        "chromsome ", seqname,
        " distal adjust done at ", date(), ".\n"
      )
    }
    chr.abun <- search_proximalCPs(
      chr.abun, curr_UTR,
      window_size, MINSIZE,
      cutEnd, search_point_START,
      search_point_END,
      filter.last
    )

    ## Step 3: search proximal CP sites
    if (!silence) {
      message(
        "chromsome ", seqname,
        " proximal search start at ", date(), ".\n"
      )
    }

    if (!silence) {
      message(
        "chromsome ", seqname,
        " proximal searched  done at ", date(), ".\n"
      )
    }

    ## Step 4: adjust proximal CP sites
    if (!silence) {
      message(
        "chromsome ", seqname,
        " proximal adjust start at ", date(), ".\n"
      )
    }
    chr.abun <- adjust_proximalCPs(
      chr.abun, PolyA_PWM,
      genome, classifier,
      classifier_cutoff,
      shift_range,
      search_point_START,
      step
    )
    if (!silence) {
      message(
        "chromsome ", seqname,
        " proximal adjust done at ", date(), ".\n"
      )
    }  
    chr.abun
  }

  if (.Platform$OS.type == "windows") {
    chr.abun <- get_CPsites(
      seqname,
      chr.cov.merge,
      conn_next_utr3,
      curr_UTR,
      depth.weight,
      background,
      z2s,
      window_size,
      long_coverage_threshold,
      adjust_distal_polyA_end,
      classifier,
      classifier_cutoff,
      shift_range,
      genome,
      step,
      MINSIZE,
      search_point_START,
      search_point_END,
      cutEnd,
      filter.last,
      PolyA_PWM,
      outdir,
      silence
    )
    ## save intermediate data before polishing
    filename <- file.path(outdir, paste0(seqname, "_CPsites.no.polishing.RDS"))
    saveRDS(chr.abun, file = filename)

    ## polishing
    chr.abun <- polish_CPs(chr.abun)
  } else {
    reg <- makeRegistry(
      file.dir = paste(outdir,
        seqname,
        sep = "_"
      ),
      conf.file = NA,
      packages = "InPAS",
      seed = 1
    )
    cluster_type <- match.arg(cluster_type)

    if (cluster_type %in% c("lsf", "sge", "slurm", "openlava", "torgue") &&
      !file.exists(template_file)) {
      stop("template_file doen't exist")
    }

    if (cluster_type == "interactive") {
      reg$cluster.functions <- makeClusterFunctionsInteractive(
        external = FALSE,
        write.logs = TRUE,
        fs.latency = 0
      )
    } else if (cluster_type == "lsf") {
      reg$cluster.functions <- makeClusterFunctionsLSF(
        template = template_file,
        scheduler.latency = 1,
        fs.latency = 65
      )
    } else if (cluster_type == "multicore") {
      reg$cluster.functions <-
        makeClusterFunctionsMulticore(
          ncpus = mc.cores,
          fs.latency = 0
        )
    } else if (cluster_type == "socket") {
      reg$cluster.functions <-
        makeClusterFunctionsSocket(
          ncpus = mc.cores,
          fs.latency = 65
        )
    } else if (cluster_type == "openlava") {
      reg$cluster.functions <- makeClusterFunctionsOpenLava(
        template = template_file,
        scheduler.latency = 1,
        fs.latency = 65
      )
    } else if (cluster_type == "sge") {
      reg$cluster.functions <- makeClusterFunctionsSGE(
        template = template_file,
        nodename = "localhost",
        scheduler.latency = 1,
        fs.latency = 65
      )
    } else if (cluster_type == "slurm") {
      reg$cluster.functions <- makeClusterFunctionsSlurm(
        template = template_file,
        array.jobs = TRUE,
        nodename = "localhost",
        scheduler.latency = 1,
        fs.latency = 65
      )
    } else if (cluster_type == "torque") {
      reg$cluster.functions <- makeClusterFunctionsTORQUE(
        template = template_file,
        scheduler.latency = 1,
        fs.latency = 65
      )
    }
    x <- ceiling(seq_along(curr_UTR) / future.chunk.size)
    curr_UTR <- split(curr_UTR, x)
    curr_UTR <- as.list(curr_UTR)
    conn_next_utr3 <- split(conn_next_utr3, x)
    chr.cov.merge <- split(chr.cov.merge, x)

    batchMap(
      fun = get_CPsites,
      args = list(
        curr_UTR = curr_UTR,
        conn_next_utr3 = conn_next_utr3,
        chr.cov.merge = chr.cov.merge
      ),
      more.args = list(
        seqname = seqname,
        depth.weight = depth.weight,
        background = background,
        z2s = z2s,
        window_size = window_size,
        long_coverage_threshold =
          long_coverage_threshold,
        adjust_distal_polyA_end =
          adjust_distal_polyA_end,
        classifier = classifier,
        classifier_cutoff = classifier_cutoff,
        shift_range = shift_range,
        genome = genome,
        step = step,
        MINSIZE = MINSIZE,
        search_point_START = search_point_START,
        search_point_END = search_point_END,
        cutEnd = cutEnd,
        filter.last = filter.last,
        PolyA_PWM = PolyA_PWM,
        outdir = outdir,
        silence = silence
      )
    )

    ## start job
    submitJobs(resources = resources)
    while (!waitForJobs(sleep = 120, timeout = Inf)) {
      Sys.sleep(120)
    }
    if (waitForJobs(sleep = 120, timeout = Inf)) {
      chr.abun <- reduceResultsList()
    }
    ## save intermediate data before polishing
    filename <- file.path(outdir, paste0(seqname, "_CPsites.no.polishing.RDS"))
    saveRDS(chr.abun, file = filename)

    ## polish
    chr.abun <- do.call("rbind", lapply(chr.abun, polish_CPs))
  }

  if (!is.null(chr.abun)) {
    utr3.shorten.UTR <- chr.utr3 %>%
      plyranges::filter(feature == "utr3") %>%
      plyranges::select(-c(feature)) %>%
      plyranges::filter(transcript %in% rownames(chr.abun))
    chr.abun <- chr.abun[utr3.shorten.UTR$transcript, , drop = FALSE]

    utr3.shorten.UTR$fit_value <- unlist(chr.abun[, "fit_value"])
    utr3.shorten.UTR$Predicted_Proximal_APA <-
      unlist(chr.abun[, "Predicted_Proximal_APA"])
    utr3.shorten.UTR$Predicted_Distal_APA <-
      unlist(chr.abun[, "Predicted_Distal_APA"])
    
    utr3.shorten.UTR$Predicted_Distal_APA_type <-
      unlist(chr.abun[, "distalCPtype"])
    utr3.shorten.UTR$Predicted_Proximal_APA_type <-
      unlist(chr.abun[, "Predicted_Proximal_APA_Type"])
    
    # utr3.shorten.UTR$type <- unlist(chr.abun[, "type"])
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

    tryCatch(
      {
        file_lock <- flock::lock(lock_filename)
        db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_db)
        res <- dbSendStatement(
          db_conn,
          paste0("DELETE FROM CPsites WHERE chr = '", seqname, "';")
        )
        dbClearResult(res)
        res <- dbSendStatement(
          db_conn,
          paste0("INSERT INTO CPsites (chr, cpsites_file)
                  VALUES ('", seqname, "','", filename, "');")
        )
        dbClearResult(res)
      },
      error = function(e) {
        print(paste(conditionMessage(e)))
      },
      finally = {
        dbDisconnect(db_conn)
        unlock(file_lock)
      }
    )
  }
  utr3.shorten.UTR
}
