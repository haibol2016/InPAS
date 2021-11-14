#' do test for dPDUI
#'
#' do test for dPDUI
#' 
#' @param eset An object of [UTR3eSet-class]. It is an output of 
#'   [get_UTR3eSet()]
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
#' @param robust logical, should the estimation of the empirical Bayes prior
#'   parameters be robustified against outlier sample variances?
#' @param ... other arguments are passed to lmFit
#' @param sqlite_db A path to the SQLite database for InPAS, i.e. the output of
#'   setup_sqlitedb().
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
#' tags <- colnames(eset@@PDUI)
#' g <- factor(gsub("\\..*$", "", tags))
#' design <- model.matrix(~ -1 + g)
#' colnames(design) <- c("Brain", "UHR")
#' contrast.matrix <- makeContrasts(contrasts = "Brain-UHR", 
#'                                  levels = design)
#' res <- test_dPDUI(eset = eset,
#'                   method = "limma",
#'                   normalize = "none",
#'                   design = design,
#'                   contrast.matrix = contrast.matrix)

test_dPDUI <- function(eset,
                      method = c("limma", "fisher.exact",
                                 "singleSample", "singleGroup"),
                      normalize = c("none", "quantiles", 
                                    "quantiles.robust",
                                    "mean", "median"),
                      design, 
                      contrast.matrix,
                      coef = 1, 
                      robust = FALSE, 
                      ...,
                      sqlite_db) {
  if (missing(eset)){
    stop("eset is required")
  }
  if (!is(eset, "UTR3eSet")){
    stop("eset is not an object of UTR3eSet class!")
  }
  
  method <- match.arg(arg = method, choices = c("limma", "fisher.exact",
                                                "singleSample", "singleGroup"))
  normalize <- match.arg(arg = normalize, choices = c("none", "quantiles", 
                                                     "quantiles.robust",
                                                     "mean", "median"))
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
      if (missing(sqlite_db)) {
        stop("sqlite_db is required for Fisher exact test.")
      }
      db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname = sqlite_db)
      metadata <- dbReadTable(db_conn, "metadata")
      dbDisconnect(db_conn)
      groups <- split(metadata$tag, metadata$condition) 
      
      if (length(groups) != 2){
        stop("The number of groups is not 2. Fisher exact",
             "test is not appropriate!")
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
  return(eset)
}
