#' use limma to analyze the PDUI
#'
#' use limma to analyze the PDUI
#'
#' @param UTR3eset an [UTR3eSet-class] object
#' @param design the design matrix of the experiment, with rows corresponding to
#'   arrays and columns to coefficients to be estimated. Defaults to the unit
#'   vector meaning that the arrays are treated as replicates. see
#'   [stats::model.matrix()]
#' @param contrast.matrix numeric matrix with rows corresponding to coefficients
#'   in fit and columns containing contrasts. May be a vector if there is only
#'   one contrast. see [limma::makeContrasts()]
#' @param coef column number or column name specifying which coefficient or
#'   contrast of the linear model is of interest. see more [limma::topTable()].
#'   default value: 1
#' @param robust logical, should the estimation of the empirical Bayes prior
#'   parameters be robustified against outlier sample variances?
#' @param ... other arguments are passed to lmFit
#'
#' @return fit results of eBayes by limma. It is an object of class
#'   [limma::MArrayLM-class] containing everything found in fit. see
#'   [limma::eBayes()]
#' @export
#' @import limma
#'
#' @seealso [singleSampleAnalysis()], [singleGroupAnalysis()], [fisherExactTest()]
#'
#' @examples
#' library(limma)
#' path <- system.file("extdata", package = "InPAS")
#' load(file.path(path, "eset.MAQC.rda"))
#' tags <- colnames(eset$PDUI.log2)
#' g <- factor(gsub("\\..*$", "", tags))
#' design <- model.matrix(~ -1 + g)
#' colnames(design) <- c("Brain", "UHR")
#' contrast.matrix <- makeContrasts(contrasts = "Brain-UHR", levels = design)
#' res <- limmaAnalysis(eset, design, contrast.matrix)
#' head(res)
limmaAnalysis <- function(UTR3eset, design, contrast.matrix,
                          coef = 1, robust = FALSE, ...) {
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
