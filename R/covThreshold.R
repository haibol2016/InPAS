#' Helper function to get intronic regions
#'
#' @param TxDb An object of [GenomicFeatures::TxDb-class]
#'
#' @return GRanges of intronic regions not overlapping exons
#' @import GenomicRanges  GenomicFeatures
#' @keywords internal

intronRegion <- function(TxDb) {
  if (!is(TxDb, "TxDb")) stop("TxDb must be an object of TxDb")
  exons <- exons(TxDb, columns = c("exon_id", "gene_id")) %>%
    plyranges::reduce_ranges_directed()
  genes <- genes(TxDb) %>%
    plyranges::reduce_ranges_directed()
  dis <- disjoin(c(genes, exons))
  ol <- findOverlaps(dis, exons)
  intronRegion <- dis[-unique(queryHits(ol))]
}


#' calculate the cutoff threshold of coverage
#'
#' calculate the cutoff threshold of coverage for long and short isoforms
#'
#' automatically determining the long_coverage_threshold by non_zero
#' intergenicRegion coverageauto determine the short_coverage_threshold by
#' quantile of non_zero intragenicRegion coverage
#' @param coverage Coverage for each sample, output of [coverageFromBedGraph()]
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param TxDb An object of [GenomicFeatures::TxDb-class]
#' @param utr3 Output of [utr3Annotation()]
#' @param chr Chromosome to be used for calculation, default is "chr1"
#' @param hugeData Is this dataset consume too much memory? if it is TRUE, the
#'   coverage will be saved into tempfiles.
#' @param groupList a list of grouped sample tag names, with the group names as 
#'   the list's name, such as list(groupA = c("sample_1", "sample_2", 
#'   "sample_3"), groupB = c("sample_4", "sample_5", "sample_6"))
#'
#' @return A numeric vector
#' @import GenomicRanges
#' @importFrom plyranges as_granges complement_ranges disjoin_ranges filter
#'   group_by mutate reduce_ranges reduce_ranges_directed remove_names select
#'   set_genome_info shift_downstream summarise
#' @keywords internal


covThreshold <- function(coverage,
                         genome,
                         TxDb,
                         removeScaffolds = FALSE,
                         utr3,
                         chr = "chr1",
                         hugeData,
                         groupList) {
  if (!is(TxDb, "TxDb")) stop("TxDb must be an object of TxDb")

  if (removeScaffolds) {
    TxDb_seqlevels <- seqlevels(TxDb)
    seqlevels(TxDb) <-
      TxDb_seqlevels[grepl("^(chr)?(\\d+|[XY])$", TxDb_seqlevels)]
  }

  totalCov <- totalCoverage(
    coverage, genome,
    hugeData, groupList
  )
  chr1totCov <- lapply(totalCov, "[[", chr)
  N <- length(chr1totCov)
  if (N > 1) {
    for (i in 2:N) {
      chr1totCov[[1]] <- chr1totCov[[1]] + chr1totCov[[i]]
    }
  }
  chr1totCov <- chr1totCov[[1]]
  intronRegion <- intronRegion(TxDb) %>%
    plyranges::filter(seqnames == chr)

  utr3Chr1Region <- utr3 %>%
    plyranges::filter(seqnames == chr & feature == "utr3")

  covBg <- function(.cvg, start, end) {
    view <- Views(.cvg, start, end)
    view <- viewApply(view, function(.ele) as.integer(.ele))
    view <- unlist(view)
    view <- view[view >= N]
    floor(quantile(view)[2])
  }
  long_coverage_threshold <- covBg(
    chr1totCov,
    start(intronRegion),
    end(intronRegion)
  )
  short_coverage_threshold <- covBg(
    chr1totCov,
    start(utr3Chr1Region),
    end(utr3Chr1Region)
  )
  ceiling(c(long_coverage_threshold / N, short_coverage_threshold / N))
}
