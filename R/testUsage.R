#' do test for dPDUI
#'
#' do test for dPDUI
#' @param CPsites outputs of [CPsites()]
#' @param coverage coverage for each sample, outputs of [coverageFromBedGraph()]
#' @param genome an object of [BSgenome::BSgenome-class]
#' @param utr3 output of [utr3Annotation()]
#' @param UTR3CDS.cov an object of [GenomicRanges::GRanges-class], output of
#'   [integrate3UTRUsage()]
#' @param hugeData logical(1), indication whether the data is huge   
#' @param BPPARAM an optional [BiocParallel::BiocParallelParam-class] instance
#'   determining the parallel back-end to be used during evaluation, or a list
#'   of BiocParallelParam instances, to be applied in sequence for nested calls
#'   to bplapply.
#' @param method the test method. see [singleSampleAnalysis()],
#'   [singleGroupAnalysis()],[fisherExactTest()], [limmaAnalysis()]
#' @param normalize the normalization method
#' @param design a design matrix of the experiment, with rows corresponding to
#'   arrays and columns to coefficients to be estimated. Defaults to the unit
#'   vector meaning that the arrays are treated as replicates. see
#'   [stats::model.matrix()]. Required for limma-based analysis.
#' @param contrast.matrix a numeric matrix with rows corresponding to
#'   coefficients in fit and columns containing contrasts. May be a vector if
#'   there is only one contrast. see [limma::makeContrasts()]. Required for 
#'   limma-based analysis.
#' @param coef column number or column name specifying which coefficient or
#'   contrast of the linear model is of interest. see more [limma::topTable()].
#'   default value: 1
#' @param robust logical, should the estimation of the empirical Bayes prior
#'   parameters be robustified against outlier sample variances?
#' @param ... other arguments are passed to lmFit
#' @param gp1 tag names involved in group 1. required for Fisher Exact test
#' @param gp2 tag names involved in group 2. required for Fisher Exact test
#'
#' @return a list with test results in a matrix.
#' @details if method is "limma", design matrix and contrast is required. if
#'   method is "fisher.exact", gp1 and gp2 is required.
#' @export
#' @seealso [singleSampleAnalysis()], [singleGroupAnalysis()],
#'   [fisherExactTest()], [limmaAnalysis()]
#' @examples
#' library(limma)
#' path <- system.file("extdata", package = "InPAS")
#' load(file.path(path, "CPs.MAQC.rda"))
#' load(file.path(path, "coverage.MAQC.rda"))
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' data(utr3.hg19)
#' tags <- names(coverage)
#' g <- factor(gsub("\\..*$", "", tags))
#' design <- model.matrix(~ -1 + g)
#' colnames(design) <- c("Brain", "UHR")
#' contrast.matrix <- makeContrasts(contrasts = "Brain-UHR", levels = design)
#' res <- testUsage(
#'   CPsites = CPs,
#'   coverage = coverage,
#'   genome = BSgenome.Hsapiens.UCSC.hg19,
#'   utr3 = utr3.hg19,
#'   method = "limma",
#'   design = design,
#'   contrast.matrix = contrast.matrix
#' )
testUsage <- function(CPsites, coverage, genome,
                      utr3, 
                      UTR3CDS.cov = NULL, 
                      hugeData = FALSE,
                      BPPARAM = NULL,
                      method = c(
                        "limma",
                        "fisher.exact",
                        "singleSample",
                        "singleGroup"
                      ),
                      normalize = c(
                        "none", "quantiles",
                        "quantiles.robust",
                        "mean", "median"
                      ),
                      design, contrast.matrix,
                      coef = 1, robust = FALSE, ...,
                      gp1, gp2) {
  stopifnot(length(CPsites) > 0 && is(CPsites, "GRanges"))
  stopifnot(length(coverage) > 0)
  if (missing(utr3) || missing(genome)) {
    stop("utr3 and genome are required.")
  }
  if (!is(genome, "BSgenome")) {
    stop("genome must be an object of BSgenome.")
  }
  if (!is(utr3, "GRanges") ||
      !all(utr3$feature %in% c("utr3", "next.exon.gap", "CDS"))) {
    stop("utr3 must be output of function of utr3Annotation")
  }
  if (hugeData && is.null(UTR3CDS.cov)){
    stop("UTR3CDS.cov, an object of GRanges containing 3' UTR coverage data for ", 
          "CPsites should be prepared before hand if it is huge data")
  }
  if (!is.null(UTR3CDS.cov) && !is(UTR3CDS.cov, "GRanges")){
    stop("UTR3CDS.cov must be an object of GRanges containing 3' UTR ", "
         coverage data for CPsites should be prepared before hand")
  }
  
  method <- match.arg(method)
  normalize <- match.arg(normalize)
  res <- NULL
  eset <- list()
  if (method != "singleSample") {
    if (!hugeData) {
      eset <- getUTR3eSet(CPsites, coverage, genome,
                          utr3, normalize,
                          BPPARAM = BPPARAM)
    } else {
      eset <- getUTR3eSetFromHugeData(CPsites, coverage, genome,
                                      utr3, UTR3CDS.cov,
                                      normalize,
                                      BPPARAM = BPPARAM)
    }

    if (method == "limma") {
      if (missing(design) || missing(contrast.matrix)) {
        stop("design and contrast.matrix are required.")
      }
      res <- limmaAnalysis(eset, design, contrast.matrix,
        coef = coef, robust = FALSE, ...
      )
    }
    if (method == "fisher.exact") {
      if (missing(gp1) || missing(gp2)) {
        stop("gp1 and gp2 are required.")
      }
      res <- fisherExactTest(eset, gp1 = gp1, gp2 = gp2)
    }
    if (method == "singleGroup") {
      res <- singleGroupAnalysis(eset)
    }
  } else {
    eset <- getUTR3eSet(CPsites, coverage, genome,
      utr3,
      BPPARAM = BPPARAM,
      singleSample = TRUE
    )
    res <- singleSampleAnalysis(eset)
  }
  eset$testRes <- res
  return(eset)
}
