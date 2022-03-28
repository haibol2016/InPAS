#' calculate the CP score
#'
#' calculate CP score by cleanUpdTSeq
#'
#' @param seqname a character(1) vector, the chromosome/scaffold's name
#' @param pos genomic positions
#' @param str DNA strand
#' @param idx offset position
#' @param idx.gp group number of the offset position
#' @param genome an object of [BSgenome::BSgenome-class]
#' @param classifier An R object for Naive Bayes classifier model, like the one
#'   in the cleanUpdTSeq package.
#' @param classifier_cutoff A numeric(1) vector. A cutoff of probability that a
#'   site is classified as true CP sites. The value should be between 0.5 and 1.
#'   Default, 0.8.
#' @return a data frame or NULL
#' @seealso [get_PAscore()]
#' @importFrom cleanUpdTSeq buildFeatureVector predictTestSet
#' @importFrom future plan multicore multisession
#' @importFrom future.apply future_lapply
#' @import GenomicRanges
#' @keywords internal
#' @author Jianhong Ou, Haibo Liu

get_PAscore2 <- function(seqname,
                         pos,
                         str,
                         idx,
                         idx.gp,
                         genome,
                         classifier,
                         classifier_cutoff) {
  if (length(pos) < 1) {
    return(NULL)
  }
  coor <- paste(seqname, pos, str, sep = "_")
  gr <- GRanges(seqname,
    IRanges(pos, pos, names = coor),
    strand = str
  )

  gr$id <- 1:length(gr)
  coor.id <- !duplicated(coor)
  gr$duplicated <- gr$id[coor.id][match(coor, coor[coor.id])]
  gr.s <- gr[coor.id]

  nbc_scoring <- function(.gr.s) {
    testSet.NaiveBayes <- buildFeatureVector(
      .gr.s,
      genome = genome,
      upstream = classifier@info@upstream,
      downstream = classifier@info@downstream,
      wordSize = classifier@info@wordSize,
      alphabet = classifier@info@alphabet,
      sampleType = "unknown",
      replaceNAdistance = 30,
      method = "NaiveBayes",
      fetchSeq = TRUE,
      return_sequences = FALSE
    )

    ## seqname is not needed
    suppressMessages({
      .pred.prob.test <-
        predictTestSet(
          testSet.NaiveBayes = testSet.NaiveBayes,
          classifier = classifier,
          outputFile = NULL,
          assignmentCutoff = classifier_cutoff,
          return_sequences = FALSE
        )
    })
    .pred.prob.test[, -2]
  }
  gr_lists <- split(gr.s, ceiling(seq_along(gr.s) / 100))
  cores <- parallelly::availableCores() -1
  scores <- future_lapply(gr_lists, nbc_scoring, 
        future.scheduling = floor(length(gr_lists)/cores))
  pred.prob.test <- do.call(rbind, scores)
  pred.prob.test <- pred.prob.test[match(
    names(gr.s),
    pred.prob.test$peak_name
  ), ,
  drop = FALSE
  ]
  if (any(duplicated(coor))) {
    ## need to recover the order of inputs
    pred.prob.test <- pred.prob.test[match(gr$duplicated, gr.s$id), ,
      drop = FALSE
    ]
    pred.prob.test[, "peak_name"] <- names(gr)
  }
  pred.prob.test <- cbind(pred.prob.test, idx, idx.gp)
  pred.prob.test <- pred.prob.test[!is.na(pred.prob.test[, "pred_class"]), ]
  pred.prob.test <- pred.prob.test[pred.prob.test[, "pred_class"] == 1, ]
  if (nrow(pred.prob.test) > 0) {
    pred.prob.test <- pred.prob.test[order(
      pred.prob.test[, "idx.gp"],
      -pred.prob.test[, "prob_true_pA"]
    ), ]
    pred.prob.test <- pred.prob.test[!duplicated(pred.prob.test[, "idx.gp"]), ]
    pred.prob.test
  } else {
    # message("no CP sites!")
    NULL
  }
}
