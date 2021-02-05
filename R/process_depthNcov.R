#' Format depth and coverage data
#' 
#' Glean sample depth data together and convert sample-wise coverage data to 
#'    chromosome-wise coverage data
#'
#' @param coverage A vector of character(n) containing paths to RDS files 
#'    storing Rle instances and sample total read depth for hugeData, or 
#'    a list containing Rle instances and sample total read depth 
#' @param outdir A path with write permission. If it is NULL, the output will
#'    be written to the same directory as the input coverage files for hugeData
#'
#' @return
#' @export
#'


process_depthNcov <- function(coverage){
    if (missing(coverage) || length(coverage) < 1){
        stop("At least one coverage is required!")
    }
    hugeData <- is.character(coverage[1])
    if (hugeData && !any(file.exists(coverage))){
        stop("Some coverage files do not exist!")
    }
    
    cov <- list()
    depth <- vector("integer", length(num_samples))
    
    if (hugeData){
        for (i in seq_along(coverage)){
            cvgNdepth <- readRDS(coverage[i])
            depth[i] <- cvgNdepth$depth
            names(depth)[i] <- names(cvgNdepth$depth)
            tag <- names(cvgNdepth)[1]
            chr_names <- names(cvgNdepth[[1]])
            for (j in chr_names){
                cov[[j]][[tag]] <- cvgNdepth[[1]][[j]]
            }
        }
        
        if (is.null(outdir)){
            outdir <- dirname(coverage[1])
        } else {
            if (!dir.exists(outdir)){
                dir.create(outdir, recursive = TRUE)
            }
        }
        outdir <- normalizePath(outdir)
        
        ## split coverage list per chromosome
        filenames <- mapply(function(.cov, .chr){
            filename <- file.path(outdir, 
                                  paste0(.chr, "_", 
                                basename(tempfile(pattern = "RleCov_"))))
            
            ## A list of Rle for a chromosome of each sample
            saveRDS(.cov, file = filename)
            filename
        }, cov, names(cov), SIMPLIFY = FALSE)
        list(cov = filenames, depth = depth)
    } else { 
        for (i in seq_along(coverage)){
            cvgNdepth <- coverage[[i]]
            depth[i] <- cvgNdepth$depth
            names(depth)[i] <- names(cvgNdepth$depth)
            tag <- names(cvgNdepth)[1]
            chr_names <- names(cvgNdepth[[1]])
            for (j in chr_names){
                cov[[j]][[tag]] <- cvgNdepth[[1]][[j]]
            }
        }
        list(cov = cov, depth = depth)
    }
}