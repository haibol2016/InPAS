#' Visualization of MSE profiles, 3' UTR coverage and minimal MSE distribution
#'
#' @param CPs A list, output from [search_proximalCPs()] or [adjust_distalCPs()] 
#'   or [adjust_proximalCPs()]
#' @param outdir A character(1) vector, specifying the output directory
#' @param MSE.plot A character(1) vector, specifying a PDF file name for 
#'   outputting plots of MSE profiles. No directory path is allowed.
#' @param coverage.plot A character(1) vector, specifying a PDF file name for
#'   outputting per-sample coverage profiles. No directory path is allowed.
#' @param min.MSE.to.end.distr.plot A character(1) vector, specifying a PDF 
#'   file name for outputting histograms showing minimal MSE distribution 
#'   relative to longer 3' UTR end. No directory path is allowed.
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom graphics hist legend
#' @importFrom stats density 
#' @export
find_minMSEDistr <- function(CPs, outdir = NULL, 
                             MSE.plot = "MSE.pdf",
                             coverage.plot = "coverage.pdf",
                             min.MSE.to.end.distr.plot = 
                                 "min.MSE.to.end.distr.pdf")
{
    if (is.null(CPs[["chr.cov.merge"]]) || 
        is.null(CPs[["fit_value"]]))
    {
        stop("CPs is not a valid list.\n", 
             "per-sample coverage data and fitted values, MSEs\n",
             "should be stored in CPs!")
    }
    if (is.null(outdir)) 
    {
        stop("outdir is needed!")
    }
    if (!dir.exists(outdir))
    {
        dir.create(outdir, recursive = TRUE)
    }
    
    files <- c(MSE.plot, coverage.plot, min.MSE.to.end.distr.plot)
    if (!all(grepl("\\.pdf$", files)) || 
        any(grepl("/", files)) ||
        !all(file.create(file.path(outdir, files))))
    {
        stop("valid .pdf file names for outputting plots of MSE profiles,\n",
             "per-sample read coverage across 3' UTRs, and distribution of\n",
             "minimal non-zero MSE to longer 3 'UTR ends are needed!")
    }
    
    fit_value <- CPs$fit_value
    chr.cov.merge <- CPs$chr.cov.merge
    flag <- !unlist(lapply(fit_value, is.null))
    
    ## plot MSE profiles
    pdf(file.path(outdir, MSE.plot), height = 5, width = 15)
    null <- mapply(function(x, y){
        plot(x, main = y,
             type = "l",
             xlab = "Relative distance to start of 3' UTR",
             ylab = "MSE")
    }, fit_value[flag], names(fit_value[flag]))
    dev.off()
    
    ## plot coverage profiles
    pdf(file.path(outdir, coverage.plot), height = 5, width =15)
    null <- mapply(function(x, y)
    {
        color <- rainbow(ncol(x))
        plot(x[, 1], main = y, 
             xlab = "Relative distance to start of 3' UTR",
             ylab = "Coverage",
             type = "l",
             col = color[1],
             ylim = c(min(pmin(x)), max(pmax(x))))
        if (ncol(x) > 1)
        {
            for (i in 2:ncol(x))
            {
                lines(1:nrow(x), x[,i], col  = color[i])
            }
        }
        legend("topright", legend=colnames(x),
               col=color, lty=1, cex=0.8,
               title="Condition", 
               text.font=4)
    }, chr.cov.merge[flag], names(chr.cov.merge[flag]))
    dev.off()
    
    ## plot minMSE distribution
    rel_dist <- unlist(lapply(fit_value[flag], function(.x) {
        length(.x) - which(.x != 0 & .x== min(.x[.x !=0]))
    }))
    
    pdf(file.path(outdir, min.MSE.to.end.distr.plot), 
        height = 6, width = 10)
    hist(rel_dist, freq = FALSE, 
         breaks = 100, xlim = c(0, 2500),
         xlab = "Relative distance to end of 3' UTR",
         ylab = "Density")
    lines(density(rel_dist), xlim = c(0, 2500),
          col = "red")
    dev.off()
}
