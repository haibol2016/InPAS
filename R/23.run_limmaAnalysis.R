#' use limma to analyze the PDUI
#'
#' use limma to analyze the PDUI
#'
#' @param UTR3eset An object of [UTR3eSet-class], output of [get_UTR3eSet()]
#' @param design A design matrix of the experiment, with rows corresponding to
#'   arrays and columns to coefficients to be estimated. Defaults to the unit
#'   vector meaning that the arrays are treated as replicates. see
#'   [stats::model.matrix()]
#' @param contrast.matrix A numeric matrix with rows corresponding to coefficients
#'   in fit and columns containing contrasts. May be a vector if there is only
#'   one contrast. see [limma::makeContrasts()]
#' @param coef An integer(1) vector specifying which coefficient or a character(1)
#'   vector specifying which contrast of the linear model is to test. see more 
#'   [limma::topTable()]. Default, 1.
#' @param robust A logical(1) vector,indicating whether the estimation of the
#'   empirical Bayes prior parameters be robustified against outlier sample
#'   variances?
#' @param ... other arguments which are passed to [limma::lmFit()]
#'
#' @return fit results of eBayes by limma. It is an object of class
#'   [limma::MArrayLM-class] containing everything found by fit. see
#'   [limma::eBayes()]
#' @export
#' @import limma
#' @author Jianhong Ou
#' @seealso [run_singleSampleAnalysis()], [run_singleGroupAnalysis()], [run_fisherExactTest()]
#'
#' @keywords internal

run_limmaAnalysis <- function(UTR3eset, 
                          design, 
                          contrast.matrix,
                          coef = 1, 
                          robust = FALSE,
                          ...) {
  if (missing(design) || missing(contrast.matrix)) {
    stop("desing and contrast.matrix is required.")
  }
  if (!is(design, "matrix")) {
    stop("design must be an design matrix")
  }
  if (!is(contrast.matrix, "matrix")) {
    stop("contrast.matrix must be an object of matrix")
  }
  if (!is(UTR3eset, "UTR3eSet")) {
    stop("UTR3eset must be an object of UTR3eSet")
  }
  short <- UTR3eset@short
  long <- UTR3eset@long
  if (!identical(rownames(short), rownames(long))) {
    stop("the rownames are not identical for short form and long form")
  }
  txid <- rep(rownames(short), 2)
  formid <- c(
    paste(rownames(short), "short", sep = ":"),
    paste(rownames(long), "long", sep = ":")
  )
  y <- rbind(short, long)
  rownames(y) <- formid
  lib.size <- colSums(y)
  y <- log2(y + 0.001)
  y <- normalizeBetweenArrays(y)
  fit <- lmFit(y, design, ...)
  fit <- contrasts.fit(fit, contrast.matrix)

  ex <- diffSplice(fit,
    geneid = txid, exonid = formid,
    robust = robust, verbose = FALSE
  )
  ts <- topSplice(ex, coef = coef, test = "simes", number = nrow(ex))
  p <- ts[match(rownames(short), ts$GeneID), "P.Value"]
  BH <- ts[match(rownames(short), ts$GeneID), "FDR"]
  out <- cbind(P.Value = p, adj.P.Val = BH)
  rownames(out) <- rownames(short)
  out
}
