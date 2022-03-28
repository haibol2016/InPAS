checkCPs <- function(CP, proximal, distal, type) {
  if (!missing(proximal)) {
    checkEquals(CP[1]$Predicted_Proximal_APA, proximal)
  }
  if (!missing(distal)) {
    checkEquals(CP[1]$Predicted_Distal_APA, distal)
  }
  if (!missing(type)) {
    checkEquals(CP[1]$Predicted_Distal_APA_type, type)
  }
}

generateGRs <- function(starts, ends, scores) {
  GRanges("chr1",
    IRanges(start = starts, end = ends),
    strand = "*", score = scores
  )
}

generateUTR3 <- function(normal = TRUE) {
  if (normal) {
    chr1_gr <- GRanges("chr1", IRanges(c(100, 600), c(599, 2000),
      names = c(
        "transcript1|gene1|utr3",
        "transcript1|gene1|next.exon.gap"
      )
    ),
    strand = "+",
    feature = c("utr3", "next.exon.gap"),
    annotatedProximalCP = "unknown",
    exon = "transcript1", 
    transcript = "transcript1|",
    gene = "1", 
    truncated = FALSE,
    symbol = "gene1"
    )
    chr1_gr
  } else {
    data(utr3.mm10)
    utr3.mm10 <- utr3.mm10[seqnames(utr3.mm10) == "chr1" & end(utr3.mm10) < 10000000]
    utr3.mm10
  }
}

genome <- BSgenome.Mmusculus.UCSC.mm10
TxDb <- TxDb.Mmusculus.UCSC.mm10.knownGene
addLockName(filename = tempfile())
outdir <- tempdir()
addInPASOutputDirectory(outdir)

get_CPs <- function(tags, conditions, filenames,
                    genome, TxDb, utr3) {
  metadata <- data.frame(
    tag = tags,
    condition = conditions,
    bedgraph_file = filenames
  )
  write.table(metadata,
    file = file.path(outdir, "metadata.txt"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  sqlite_db <- setup_sqlitedb(
    metadata = file.path(
      outdir,
      "metadata.txt"
    ),
    outdir
  )
  cov <- list()
  for (i in seq_along(filenames)) {
    cov[[tags[i]]] <- get_ssRleCov(
      bedgraph = filenames[i],
      tag = tags[i],
      genome = genome,
      sqlite_db = sqlite_db,
      outdir = outdir,
      chr2exclude = NULL
    )
  }

  data4CPsitesSearch <- setup_CPsSearch(sqlite_db,
    genome = genome,
    chr.utr3 = utr3,
    seqname = "chr1",
    background = "same_as_long_coverage_threshold",
    TxDb = TxDb,
    hugeData = TRUE,
    outdir = outdir,
    coverage_threshold = 5,
    silence = FALSE
  )

  CP <- search_CPs(
    seqname = "chr1",
    sqlite_db = sqlite_db,
    genome = genome,
    MINSIZE = 10,
    window_size = 100,
    search_point_START = 50,
    search_point_END = NA,
    cutEnd = 0,
    filter.last = FALSE,
    adjust_distal_polyA_end = TRUE,
    long_coverage_threshold = 2,
    PolyA_PWM = NA,
    classifier = NA,
    classifier_cutoff = 0.8,
    shift_range = 100,
    step = 1,
    outdir = outdir,
    silence = FALSE
  )

  list(CP = CP, sqlite_db = sqlite_db)
}

NOtest_searchDistalCPs <- function() {
  # drop from end
  utr3 <- generateUTR3()
  testGR <- generateGRs(
    c(5, 600, 1001),
    c(599, 1000, 1500), c(40, 10, 1)
  )

  filename <- tempfile()
  export(testGR, filename, format = "BEDGraph")
  CP <- get_CPs(
    tags = "test", conditions = "test",
    filenames = filename, genome = genome,
    TxDb = TxDb, utr3 = utr3
  )$CP
  checkCPs(CP, 599, 1000, "extended novel distal")

  # local background
  set.seed(393993)
  utr3 <- generateUTR3(normal = FALSE)
  background <- rep(rnbinom(140000, 10e9, mu = 2), each = 50)
  signal <- rep(
    c(0, 5, 10, 11, 0),
    c(1909199, 276, 998, 189, 5089338)
  ) + background
  signal <- signal[1899200:1920662]
  signal <- rle(signal)
  signal.cumsum <- cumsum(signal$lengths)
  testGR <- generateGRs(
    4899200 + c(0, signal.cumsum[-length(signal.cumsum)]),
    4899199 + signal.cumsum, signal$values
  )
  export(testGR, filename, format = "BEDGraph")

  CP <- get_CPs(
    tags = "test", conditions = "test",
    filenames = filename, genome = genome,
    TxDb = TxDb, utr3 = utr3
  )$CP

  checkCPs(CP, proximal = 4909451, type = "extended novel distal")
}

test_CPsites_utr3Usage <- function() {
  utr3 <- generateUTR3()
  testGR <- generateGRs(c(5, 600),
    c(599, 1000),
    scores = c(40, 10)
  )

  filename <- tempfile()
  export(testGR, filename, format = "BEDGraph")

  CP <- get_CPs(
    tags = "test", conditions = "test",
    filenames = filename, genome = genome,
    TxDb = TxDb, utr3 = utr3
  )

  checkCPs(CP$CP, 599, 1000, "extended novel distal")

  sqlite_db <- CP$sqlite_db
  utr3_cds_cov <- get_regionCov(
    chr.utr3 = utr3,
    sqlite_db,
    outdir,
    phmm = FALSE
  )

  res <- get_UTR3eSet(sqlite_db,
    normalize = "none",
    singleSample = FALSE
  )

  testGR2 <- GRanges("chr1", IRanges(5, 1000),
    strand = "*", score = 20
  )

  filename2 <- tempfile()
  export(testGR2, filename2, format = "BEDGraph")

  CP <- get_CPs(
    tags = c("test1", "test2"),
    conditions = c("test1", "test2"),
    filenames = c(filename, filename2),
    genome = genome, TxDb = TxDb, utr3 = utr3
  )
  checkCPs(CP$CP, 599, 1000, "extended novel distal")

  sqlite_db <- CP$sqlite_db
  utr3_cds_cov <- get_regionCov(
    chr.utr3 = utr3,
    sqlite_db,
    outdir,
    phmm = FALSE
  )

  res <- get_UTR3eSet(sqlite_db,
    normalize = "none",
    singleSample = FALSE
  )
  checkEqualsNumeric(res@PDUI[, "test1"], 0.25, tolerance = 1.0e-2)
  checkEqualsNumeric(res@PDUI[, "test2"], 1, tolerance = 1.0e-2)

  test_out <- InPAS:::run_fisherExactTest(
    UTR3eset = res,
    gp1 = "test1",
    gp2 = "test2"
  )
  res$testRes <- test_out
  res <- filter_testOut(res = res, gp1 = "test1", gp2 = "test2")
  checkEquals(as.logical(res$PASS), TRUE)
}

test_getCov_seqlevelsStyle <- function() {
  utr3 <- generateUTR3()
  testGR <- generateGRs(
    c(5, 400, 1001),
    c(399, 1000, 1500),
    c(40, 10, 1)
  )

  seqlevelsStyle(testGR) <- "NCBI"

  filename <- tempfile()
  export(testGR, filename, format = "BEDGraph")

  err <- try({
    CP <- get_CPs(
      tags = c("test1"),
      conditions = c("test1"),
      filenames = c(filename),
      genome = genome, TxDb = TxDb,
      utr3 = utr3
    )
  })
  checkEquals(is(err, "try-error"), TRUE)
}
