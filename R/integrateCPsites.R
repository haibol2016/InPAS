#' Combine and format CPsites
#'
#' @param shorten_UTR_estimation_files A vector of character, containing paths to all
#'   the outputs of [estimateCPsites()] of an experiment. absolute paths or
#'   paths relative to current work directory
#' @param utr3 An object of [GenomicRanges::GRanges-class], output from the 
#'   [utr3Annotation()]
#'
#' @return An object of \code{\link[GenomicRanges]{GRanges}}
#' @export

integrateCPsites <- function(shorten_UTR_estimation_files, 
                             utr3){
    if (missing(shorten_UTR_estimation_files)){
        stop("shorten_UTR_estimation_files is required")
    }
    if(!any(file.exists(shorten_UTR_estimation_files))){
        stop("not all shorten_UTR_estimation_files exist")
    }
    if (missing(utr3)){
        stop("utr3 is required")
    }
    if (!is(utr3, "GRanges")){
        stop("utr3 must be a GRanges object")
    }
    len <- length(shorten_UTR_estimation_files)
    shorten_UTR_estimation <- vector("list", len)
    
    for (i in seq_along(shorten_UTR_estimation_files)){
        shorten_UTR_estimation[[i]] <- readRDS(shorten_UTR_estimation_files[i])
    }
    shorten_UTR_estimation <-
        do.call(
            rbind,
            unname(shorten_UTR_estimation[!sapply(
                shorten_UTR_estimation,
                is.null
            )])
        )
    utr3.shorten.UTR <- utr3 %>%
        plyranges::filter(feature == "utr3") %>%
        plyranges::select(-c(feature)) %>%
        plyranges::filter(transcript %in%
                              rownames(shorten_UTR_estimation))
    
    shorten_UTR_estimation <-
        shorten_UTR_estimation[utr3.shorten.UTR$transcript, , drop = FALSE]
    
    utr3.shorten.UTR$fit_value <- unlist(shorten_UTR_estimation[, "fit_value"])
    utr3.shorten.UTR$Predicted_Proximal_APA <-
        unlist(shorten_UTR_estimation[, "Predicted_Proximal_APA"])
    utr3.shorten.UTR$Predicted_Distal_APA <-
        unlist(shorten_UTR_estimation[, "Predicted_Distal_APA"])
    utr3.shorten.UTR$type <- unlist(shorten_UTR_estimation[, "type"])
    utr3.shorten.UTR$utr3start <- unlist(shorten_UTR_estimation[, "utr3start"])
    utr3.shorten.UTR$utr3end <- unlist(shorten_UTR_estimation[, "utr3end"])
    utr3.shorten.UTR <-
        utr3.shorten.UTR[!is.na(utr3.shorten.UTR$Predicted_Proximal_APA)]
}
