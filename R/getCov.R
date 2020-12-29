#' extract coverage from bedgraph file
#'
#' extract coverage from bedgraph file
#'
#' @param bedgraph bedGraph file names
#' @param genome an object [BSgenome::BSgenome-class]
#' @param seqLen a named vector containing lengths of each chromosome. output
#'   from [seqLen()]
#' @param removeScaffolds logical(1), whether the scaffolds should be removed
#'   from the genome
#'
#' @return an object of [S4Vectors::Rle-class] for a sample coverage
#' @importFrom dplyr as_tibble mutate filter arrange bind_rows group_by
#'   left_join summarise n bind_cols syms desc pull
#' @importFrom GenomeInfoDb seqlengths seqlevelsStyle mapSeqlevels
#' @import readr
#' @import S4Vectors GenomicRanges
#' @importFrom BSgenome getSeq matchPWM
#' @importFrom IRanges IRanges quantile viewApply viewMeans quantile
#' @importFrom magrittr %>%
#' @seealso [coverageFromBedGraph()]


getCov <- function(bedgraph, genome,
                   seqLen,
                   removeScaffolds = FALSE) {
  seqnames.bedfile <-
    as.data.frame(read_tsv(
      bedgraph,
      comment = "#",
      col_names = FALSE, skip = 0,
      col_types = cols(
        X1 = col_factor(),
        X2 = col_skip(),
        X3 = col_skip(),
        X4 = col_skip()
      )
    ))[, 1]

  ## filter out scaffolds and mitochondrial genome if specified
  seqnames <- trimSeqnames(genome, removeScaffolds)

  seqStyle <- seqlevelsStyle(genome)
  seqStyle.bed <- seqlevelsStyle(levels(seqnames.bedfile))

  ## if seqlevel style different between BSgenome and bedgraph file,
  ## make them the same. This is not safe-proof.

  if (!any(seqStyle.bed == seqStyle)) {
    message("seqlevelsStyle of genome is different from bedgraph file.")
    ## convert to seqStyled seqnames if matched, otherwise character_NA_: A names vector
    levels <-
      mapSeqlevels(levels(seqnames.bedfile), seqStyle)
    if (!is.character(levels)) {
      id <- apply(levels, 1, function(.ele) {
        sum(seqnames %in% .ele)
      })
      levels <- levels[id == max(id), ][1, ]
    }
    levels[is.na(levels)] <- names(levels)[is.na(levels)]
    levels(seqnames.bedfile) <- levels
  }
  seqnames <- sort(intersect(levels(seqnames.bedfile), seqnames))
  if (length(seqnames) < 1) {
    stop(paste(
      "there is no intersect seqname in",
      bedgraph, "and genome"
    ))
  }

  summaryFunction <- function(seqname) {
    seqL <- seqLen[seqname]

    ## apply some trick here by using "FALSE"
    lines2read <- Rle(c(FALSE, seqnames.bedfile == seqname))
    true <- which(runValue(lines2read))
    skip <- runLength(lines2read)[true - 1]
    skip[1] <- skip[1] - 1
    nrow <- runLength(lines2read)[true]

    ## To be efficient, sort the bedgraph by chromosomes first
    ## If bedgraph files are generated using bedtools genomecov,
    ## Then this is the case
    dat <- read_tsv(bedgraph,
      comment = "#",
      col_names = FALSE,
      skip = skip[1],
      n_max = nrow[1],
      col_types = cols(
        X1 = col_skip(),
        X2 = col_integer(),
        X3 = col_integer(),
        X4 = col_integer()
      )
    )

    if (length(true) > 1) {
      cul_skip <- skip[1] + nrow[1]
      for (i in 2:length(true)) {
        cul_skip <- cul_skip + skip[i]
        lines <-
          read_tsv(bedgraph,
            comment = "#",
            col_names = FALSE,
            skip = cul_skip,
            n_max = nrow[i],
            col_types = cols(
              X1 = col_skip(),
              X2 = col_integer(),
              X3 = col_integer(),
              X4 = col_integer()
            )
          )
        dat <- dplyr::bind_rows(dat, lines)
        cul_skip <- cul_skip + nrow[i]
      }
    }

    # convert bedgraph 0-based index to 1-based index for GRanges
    dat <- dat %>% dplyr::mutate(X2 = X2 + 1)

    ## convert uncovered regions as gaps
    ## This step can be avoided when generate bedgraph with bedtools:
    ## bedtools genomecov -bga -split -ibam  <coordinate-sorted BAM file>
    gaps <-
      as_tibble(as.data.frame(gaps(IRanges(
        dplyr::pull(dat, 1),
        dplyr::pull(dat, 2)
      ),
      start = 1, end = seqL
      )))
    if (nrow(gaps) > 0) {
      gaps <- dplyr::bind_cols(gaps[, 1:2],
        V4 = rep(0, nrow(gaps))
      )
      colnames(gaps) <- colnames(dat)
      dat <- dplyr::bind_rows(dat, gaps)
    }
    dat <- dat %>%
      dplyr::arrange(X2) %>%
      dplyr::mutate(width = X3 - X2 + 1)
    # create a RLe object
    Rle(dplyr::pull(dat, 3), dplyr::pull(dat, width))
  }
  cov <- lapply(seqnames, summaryFunction)
  names(cov) <- seqnames
  cov
}
