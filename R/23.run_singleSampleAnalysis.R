#' do APA analysis for a single sample
#'
#' @description do APA event analysis for a single sample Using Poisson Hidden
#'   Markov models
#' @details the test will be performed by comparing a two-state versus an
#'   one-state Poisson Hidden Markov models.
#'
#' @param UTR3eset the output of [get_UTR3eSet()]
#' @importFrom depmixS4 depmix fit posterior getpars
#' @importFrom stats BIC fft fisher.test loess.smooth model.matrix p.adjust
#'   pnorm poisson t.test
#' @return a matrix containing test results
#' @seealso [UTR3eSet-class], [get_UTR3eSet()], [depmixS4::depmix()]
#' @keywords internal
#' @author Jianhong Ou
#' @examples
#' path <- system.file("extdata", package = "InPAS")
#' load(file.path(path, "eset.MAQC.rda"))
#' res <- InPAS:::run_singleSampleAnalysis(eset)

run_singleSampleAnalysis <- function(UTR3eset) {
  ## could we use coverage instead of counts for poisson distribution?
  phmm <- function(counts) {
    counts <- floor(counts)
    at <- counts[1]
    counts <- counts[-1]
    mu.1 <- mean(counts[1:at]) + 1
    mu.2 <- mean(counts[-(1:at)]) + 1
    m1 <- depmix(
      response = counts ~ 1, nstates = 1,
      data = as.data.frame(counts),
      family = poisson(),
      respstart = log(mean(counts))
    )
    m2 <- depmix(
      response = counts ~ 1, nstates = 2,
      data = as.data.frame(counts),
      family = poisson(),
      respstart = c(log(mu.1), log(mu.2))
    )
    fm1 <- try(fit(m1, verbose = FALSE), silent = TRUE)
    if (inherits(fm1, "try-error")) {
      return(1)
    }
    fm2 <- try(fit(m2, verbose = FALSE), silent = TRUE)
    if (inherits(fm2, "try-error")) {
      return(1)
    }
    if (BIC(fm2) > BIC(fm1)) {
      return(1)
    }
    sts <- posterior(fm2)$state
    pars <- getpars(fm2)
    hi.state <- ifelse(pars[7] > pars[8], 1, 2)
    keytrans <- ifelse(pars[7] > pars[8], pars[5], pars[4])
    if (sts[1] != hi.state) {
      return(1)
    }
    if (sum(sts == sts[1]) - at > 1) {
      return(1)
    }
    return(keytrans)
  }

  P.Value <- lapply(UTR3eset$signals, phmm)
  adj.P.Val <- p.adjust(P.Value, method = "BH")
  short <- UTR3eset$short[, 1]
  long <- UTR3eset$long[, 1]
  PDUI <- long / (long + short)
  cbind(
    short.mean = short,
    long.mean = long,
    PDUI, P.Value, adj.P.Val
  )
}
