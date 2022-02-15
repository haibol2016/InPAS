#' do test for dPDUI
#'
#' do test for dPDUI
#'
#' @param eset An object of [UTR3eSet-class]. It is an output of
#'   [get_UTR3eSet()]
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   setup_sqlitedb().
#' @param outdir A character(1) vector, a path with write permission for storing
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @param method A character(1), indicating the method for testing dPDUI. It can
#'   be "limma", "fisher.exact", "singleSample", or "singleGroup"
#' @param normalize A character(1), indicating the normalization method. It can
#'   be "none", "quantiles", "quantiles.robust", "mean", or "median"
#' @param design a design matrix of the experiment, with rows corresponding to
#'   samples and columns to coefficients to be estimated. Defaults to the unit
#'   vector meaning that the samples are treated as replicates. see
#'   [stats::model.matrix()]. Required for limma-based analysis.
#' @param contrast.matrix a numeric matrix with rows corresponding to
#'   coefficients in fit and columns containing contrasts. May be a vector if
#'   there is only one contrast. see [limma::makeContrasts()]. Required for
#'   limma-based analysis.
#' @param coef column number or column name specifying which coefficient or
#'   contrast of the linear model is of interest. see more [limma::topTable()].
#'   default value: 1
#' @param robust A logical(1) vector, indicating whether the estimation of the empirical Bayes prior
#'   parameters should be robustified against outlier sample variances.
#' @param ... other arguments are passed to lmFit
#'
#' @return An object of [UTR3eSet-class], with the last element \code{testRes} containing the test results in a matrix.
#' @details if method is "limma", design matrix and contrast is required. if
#'   method is "fisher.exact", gp1 and gp2 is required.
#' @export
#' @seealso [run_singleSampleAnalysis()], [run_singleGroupAnalysis()],
#'   [run_fisherExactTest()], [run_limmaAnalysis()]
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' library(limma)
#' path <- system.file("extdata", package = "InPAS")
#' load(file.path(path, "eset.MAQC.rda"))
#' tags <- colnames(eset@PDUI)
#' g <- factor(gsub("\\..*$", "", tags))
#' design <- model.matrix(~ -1 + g)
#' colnames(design) <- c("Brain", "UHR")
#' contrast.matrix <- makeContrasts(
#'   contrasts = "Brain-UHR",
#'   levels = design
#' )
#' res <- test_dPDUI(
#'   eset = eset,
#'   sqlite_db,
#'   method = "limma",
#'   normalize = "none",
#'   design = design,
#'   contrast.matrix = contrast.matrix
#' )
test_dPDUI <- function(eset,
                       sqlite_db,
                       outdir = getInPASOutputDirectory(),
                       method = c(
                         "limma", "fisher.exact",
                         "singleSample", "singleGroup"
                       ),
                       normalize = c(
                         "none", "quantiles",
                         "quantiles.robust",
                         "mean", "median"
                       ),
                       design,
                       contrast.matrix,
                       coef = 1,
                       robust = FALSE,
                       ...) {
  if (missing(eset)) {
    stop("eset is required")
  }
  if (!is(eset, "UTR3eSet")) {
    stop("eset must be an object of UTR3eSet")
  }
  if (!is.character(outdir) || length(outdir) != 1) {
    stop("An explicit output directory is required")
  } else {
    outdir <- file.path(outdir, "009.test_dPDUI")
    if (!dir.exists(outdir)) {
      dir.create(outdir,
        recursive = TRUE,
        showWarnings = FALSE
      )
    }
    outdir <- normalizePath(outdir)
  }
  if (!is.logical(robust) || length(robust) != 1) {
    stop("robust must be a logical(1) vector")
  }

  method <- match.arg(arg = method, choices = c(
    "limma", "fisher.exact",
    "singleSample", "singleGroup"
  ))
  normalize <- match.arg(arg = normalize, choices = c(
    "none", "quantiles",
    "quantiles.robust",
    "mean", "median"
  ))
  res <- NULL
  if (method != "singleSample") {
    if (method == "limma") {
      if (missing(design) || missing(contrast.matrix)) {
        stop("design and contrast.matrix are required.")
      }
      res <- run_limmaAnalysis(eset, design, contrast.matrix,
        coef = coef, robust = FALSE, ...
      )
    }
    if (method == "fisher.exact") {
      if (missing(sqlite_db) || length(sqlite_db) != 1 ||
        !file.exists(sqlite_db)) {
        stop("sqlite_db is required for Fisher exact test.")
      }
      db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_db)
      metadata <- dbReadTable(db_conn, "metadata")
      dbDisconnect(db_conn)
      groups <- split(metadata$tag, metadata$condition)

      if (length(groups) != 2) {
        stop(
          "The number of groups is not 2. Fisher exact",
          "test is not appropriate!"
        )
      } else {
        res <- run_fisherExactTest(eset, gp1 = groups[[1]], gp2 = groups[[2]])
      }
    }
    if (method == "singleGroup") {
      res <- run_singleGroupAnalysis(eset)
    }
  } else {
    res <- run_singleSampleAnalysis(eset)
  }
  eset$testRes <- res
  saveRDS(eset, file = file.path(outdir, "dPDUI.test.RDS"))
  return(eset)
}
