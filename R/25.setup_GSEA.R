#' prepare files for GSEA analysis
#'
#' output the log2 transformed delta PDUI txt file, chip file, rank file and
#'  phynotype label file for GSEA analysis
#'
#' @param eset A [UTR3eSet-class] object, output of [test_dPDUI()]
#' @param groupList A list of grouped sample tag names, with the group names as 
#'   the list's name, such as list(groupA = c("sample_1", "sample_2", 
#'   "sample_3"), groupB = c("sample_4", "sample_5", "sample_6"))
#' @param preranked A logical(1) vector, out preranked or not
#' @param rankBy A character(1) vector, indicating how the gene list is ranked.
#'   It can be "logFC" or "P.value".
#' @param outdir A character(1) vector, a path with write permission for storing 
#'   the files for GSEA analysis. If it doesn't exist, it will be created.
#' @param rnkFilename A character(1) vector, specifying a filename for the
#'   preranked file
#' @param chipFilename A character(1) vector, specifying a filename for the
#'   chip file
#' @param dataFilename A character(1) vector, specifying a filename for the
#'   dataset file
#' @param PhenFilename A character(1) vector, specifying a filename for the
#'   file containing samples' phenotype labels
#'   
#' @seealso  data formats for GSEA.
#' \url{https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats}
#'
#' @import Biobase
#' @importFrom utils write.table
#' @export
#' @author Jianhong Ou, Haibo Liu
#' @examples
#' if (interactive()) {
#'   file <- system.file("extdata", "eset.MAQC.rda", package = "InPAS")
#'   load(file)
#'   gp1 <- c("Brain.auto", "Brain.phiX")
#'   gp2 <- c("UHR.auto", "UHR.phiX")
#'   groupList <- list(Brain = gp1, UHR = gp2)
#'   prepare4GSEA(eset, 
#'                groupList = groupList, 
#'                outdir = tempdir(),
#'                preranked = TRUE,
#'                rankBy = "logFC")
#' }

setup_GSEA <- function(eset,
                       groupList,
                       outdir,
                       preranked = TRUE,
                       rankBy = c("logFC", "P.value"),
                       rnkFilename = "InPAS.rnk",
                       chipFilename = "InPAS.chip",
                       dataFilename = "dPDUI.txt",
                       PhenFilename = "group.cls") {
  if (!is(eset, "UTR3eSet")) {
    stop("eset must be an object of UTR3eSet class")
  }
  if (missing(outdir)){
    stop("An explicit output directory is required")
  } else {
    outdir <- file.path(outdir, "GSEA.files")
    if (!dir.exists(outdir)){
      dir.create(outdir, recursive = TRUE, 
                 showWarnings = FALSE)
    }
    outdir <- normalizePath(outdir, mustWork = TRUE)
  }
  
  if (preranked) {
    if (length(eset$testRes) < 1) {
      stop("There is no p.value for ranking. Please try test_dPDUI first")
    }
    rankBy <- match.arg(arg = rankBy, choices = c("logFC", "P.value"))
    testRes <- eset$testRes
    usage <- eset$usage
    if (!identical(usage$transcript, rownames(testRes))) {
      stop("The rownames of eset slots are not identical")
    }
    if (rankBy == "P.value"){
      # rank by p-values, why not by log2FC?
      rank <- data.frame(symbol = usage$symbol, 
                         p.value = 1 - testRes[, "P.Value"])
    } else if (rankBy == "logFC") {
      rank <- data.frame(symbol = usage$symbol, 
                         logFC = testRes[, "logFC"])
    }
    write.table(rank, file.path(outdir, rnkFilename),
      sep = "\t", quote = FALSE, row.names = F, col.names = F)
  } else {
    eset <- as(eset, "ExpressionSet")
    exprs <- exprs(eset)
    samples <- colnames(exprs)
    fD <- fData(eset)
    chip <- fD[, c("transcript", "symbol", "gene")]
    colnames(chip) <- c("Probe Set ID", "Gene Symbol", "Gene Title")
    dir.create(outdir, showWarnings = FALSE)
    write.table(chip, file.path(outdir, chipFilename),
      sep = "\t", quote = FALSE, row.names = F
    )
    exprs <- cbind(NAME = rownames(exprs), exprs)
    write.table(exprs, file.path(outdir, dataFilename),
      sep = "\t", quote = FALSE, row.names = F
    )
    if ((!missing(groupList)) && is.list(groupList)) {
      gp <- rep(names(groupList), sapply(groupList, length))
      names(gp) <- unlist(groupList)
      if (all(samples %in% names(gp))) {
        gp <- gp[samples]
        levels <- unique(gp)
        out <- c(
          paste(length(gp), length(levels), 1),
          paste("#", paste(levels, collapse = " ")),
          paste(gp, collapse = " ")
        )
        writeLines(out, con = file.path(outdir, PhenFilename))
      }
    }
  }
}
