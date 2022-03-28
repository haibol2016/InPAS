#' Adjust the proximal CP sites
#'
#' Adjust the proximal CP sites by PolyA PWM and cleanUpdTSeq. A few candidate
#' sites, which are ranked by MSE from low to high, are used as input for
#' adjusting. The final sites are the one with best score as PA sites, which are
#' not necessary from the lowest MSE sites.
#'
#' @param CPs the outputs of [search_proximalCPs()]
#' @param PolyA_PWM PolyA position weight matrix
#' @param genome a [BSgenome::BSgenome-class] object
#' @param classifier cleanUpdTSeq classifier
#' @param classifier_cutoff cutoff value of the classifier
#' @param shift_range the searching range for the better CP sites
#' @param search_point_START just in case there is no better CP sites
#' @param step An integer, specifying an adjusting step, default 1, means 
#'   adjusting by each base by cleanUpdTSeq.
#' @param DIST2ANNOAPAP An integer, specifying a cutoff for annotate MSE valleys
#'   with known proximal APAs in a given downstream distance. Default is 1500.
#' @return keep same as [search_proximalCPs()], which can be handled by
#'   [polish_CPs()].
#' @seealso [search_proximalCPs()], [polish_CPs()], [adjust_proximalCPsByPWM()],
#'   [adjust_proximalCPsByNBC()], [get_PAscore()],[get_PAscore2()]
#' @keywords internal
#' @author Jianhong Ou

adjust_proximalCPs <- function(CPs,
                               PolyA_PWM,
                               genome,
                               classifier,
                               classifier_cutoff,
                               shift_range,
                               search_point_START,
                               step = 1,
                               DIST2ANNOAPAP = 1000) {
  dCPs <- CPs$dCPs
  seqnames <- as.character(dCPs$seqnames)
  strands <- as.character(dCPs$strand)
  starts <- lapply(
      CPs$chr.cov.merge,
      function(.ele) {
          as.numeric(gsub("^.*_SEP_", "", rownames(.ele)))
      })
    
  # The line below is not correct. Don't reverse the start (By Haibo Liu).
  # See "pos <- if (strand == "+") start + idx - 1 else start - idx + 1" in "15.
  # adjust_proximalCPs.R"
  # starts[strands == "-"] <- lapply(starts[strands == "-"], rev)
  starts <- sapply(starts, `[`, 1)

  ## to be conservative, replace predicted proximal CP sites with the known CP
  ## sites within the range of 200 bp.
  saved.id <- CPs$saved.id ## list
  Predicted_Proximal_APA <-
    mapply(function(pos, known_proximal, start, strand) {
      anno.info <- list(annotated = FALSE, site = pos)
      if (!all(is.na(known_proximal))) {
        if (!all(is.na(pos))) {
          pos <- {
            if (strand == "+") {
              start + pos - 1
            } else {
              start - pos + 1
            }
          }
          annotated.sites <- vector("numeric")
          for (i in seq_along(pos)) {
            dist <- known_proximal - pos[i]
            if (strand == "+") {
              dist <- dist[dist > 0]
              if (length(dist) > 0 && any(dist <= DIST2ANNOAPAP)) {
                annotated.sites[i] <- known_proximal[dist == min(dist)][1]
              }
            } else {
              dist <- dist[dist < 0]
              if (length(dist) > 0 && any(dist >= -DIST2ANNOAPAP)) {
                annotated.sites[i] <- known_proximal[dist == max(dist)][1]
              }
            }
          }
          if (length(annotated.sites) > 0) {
            anno.info$annotated <- TRUE
            anno.info$site <- annotated.sites[1]
          } 
        } 
      }
      anno.info
    }, CPs$Predicted_Proximal_APA,
    saved.id, starts, strands,
    SIMPLIFY = FALSE
    )

  annotated_proximal_apa_flag <-
    unlist(lapply(Predicted_Proximal_APA, function(x) {
      x[[1]]
    }))

  annotated_proximal_apa_pos <- lapply(Predicted_Proximal_APA, 
    function(x) {
       x[[2]]
  })

  dCPs$Predicted_Proximal_APA_Type <- NA_character_
  dCPs$Predicted_Proximal_APA[annotated_proximal_apa_flag] <-
    annotated_proximal_apa_pos[annotated_proximal_apa_flag]
  dCPs$Predicted_Proximal_APA_Type[annotated_proximal_apa_flag] <-
    "Known proximal APA"

  ## why do adjustment for all candidate sites? not necessary
  ## should be resolved
  ## By Haibo, only adjusted not annotated proximal CPs
  CPs$Predicted_Proximal_APA[annotated_proximal_apa_flag] <- NA
  CPs$fit_value_min <- vector("list", length(CPs$fit_value))
  CPs$fit_value_min <- lapply(
      CPs$fit_value,
      function(x) {
              if (!is.null(x)) {
                x <- min(x[!is.na(x) & x != 0])
              } else {NA}
      })

  idx.list <- CPs$Predicted_Proximal_APA
  
  if (is(classifier, "PASclassifier")) {
    cov_diff.list <- CPs$fit_value
    idx.list <- InPAS:::adjust_proximalCPsByNBC(
      idx.list, cov_diff.list,
      seqnames, starts, strands,
      genome,
      classifier,
      classifier_cutoff,
      shift_range,
      search_point_START,
      step = step
    )
  } else if (is(PolyA_PWM, "matrix")) {
    ## improvement needed here because only test if one sequence around position
    ## a hit by PWM matrix with score  > 70%
    idx.list <- InPAS:::adjust_proximalCPsByPWM(
      idx.list, PolyA_PWM, seqnames, starts,
      strands, genome, shift_range,
      search_point_START
    )
  }

  ## relative coordinates of CP sites
  flag <- CPs$flag
  CPs$dCPs <- dCPs
  CPs$flag <- flag & !annotated_proximal_apa_flag
  ## if not a valid APA after adjustment, keep the pre-adjustment APA instead
  ## This is more robust.
  CPs$Predicted_Proximal_APA[flag] <- ifelse(is.na(idx.list[flag]), 
                                             CPs$Predicted_Proximal_APA[flag], 
                                             idx.list[flag])
  CPs
}
