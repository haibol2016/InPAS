#' extract long and short 3UTR region
#'
#' extract long and short 3UTR region
#'
#' @param .grs output of [search_CPs()]
#'
#' @return A [GenomicRanges::GRanges-class] object with short form and long 3' UTR forms
#' @keywords internal
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @author Jianhong Ou

get_UTR3region <- function(.grs) {
  .ra <- c(
    ranges(.grs),
    IRanges(min(
      .grs$Predicted_Proximal_APA[1],
      .grs$Predicted_Distal_APA[1]
    ),
    max(
      .grs$Predicted_Proximal_APA[1],
      .grs$Predicted_Distal_APA[1]
    ),
    names = "CPsite"
    )
  )
  .ra <- disjoin(.ra)
  .ra <- .ra[width(.ra) > 1]
  if (as.character(strand(.grs))[1] == "+") {
    short <- .ra[start(.ra) < .grs$Predicted_Proximal_APA[1], ]
    long <- .ra[start(.ra) >= .grs$Predicted_Proximal_APA[1], ]
  } else {
    long <- .ra[start(.ra) < .grs$Predicted_Proximal_APA[1], ]
    short <- .ra[start(.ra) >= .grs$Predicted_Proximal_APA[1], ]
  }
  long <- reduce(long)
  short <- reduce(short)
  if (length(long) == 0) {
    gr <- GRanges(as.character(seqnames(.grs))[1], short,
      strand = as.character(strand(.grs))[1],
      source = rep("short", length(short))
    )
  } else {
    if (length(short) == 0) {
      gr <- GRanges(as.character(seqnames(.grs))[1], long,
        strand = as.character(strand(.grs))[1],
        source = rep("long", length(long))
      )
    } else {
      gr <- c(
        GRanges(as.character(seqnames(.grs))[1], long,
          strand = as.character(strand(.grs))[1],
          source = rep("long", length(long))
        ),
        GRanges(as.character(seqnames(.grs))[1], short,
          strand = as.character(strand(.grs))[1],
          source = rep("short", length(short))
        )
      )
    }
  }
  for (mc in colnames(mcols(.grs))) {
    if (mc != "exon") {
      mcols(gr)[, mc] <- mcols(.grs)[1, mc]
    }
  }
  gr
}
