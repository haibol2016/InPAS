#' Process coverage files for hugData
#' 
#' Process per-sample coverage files of an experiment into a list
#'   per chromosome coverage file and depth per sample
#'
#' @param coverageFiles A vector of character(n), paths to R objects storing
#'   coverage filenames and depth for each sample, which are output of
#'   bedgraph2RleCov() of Bam2RleCov() with hugeData set to "TRUE"
#' @param tags unique sample name tags for each corresponding bedGraph file
#'
#' @return A list of coverage with tags as names
#' @export


getCov4HugeData <- function(coverageFiles, outdir){
    if (missing(coverageFiles) || length(coverageFiles) < 1) {
        stop("coverageFiles of length >= 1 are required!")
    }
    if (!any(file.exists(coverageFiles))) {
        stop("Not all coverageFiles exist!")
    }
    if (missing(outdir)){
        stop("A directory of write permission is required!")
    }
    
    od <- file.path(outdir, "chromosme-wise.RleCov")
    if (!dir.exists(od)){
        dir.create(od, recursive = TRUE)
    }
    
    ## extract sample read depth and rearrange sample-wise coverage filenames 
    ## into chromosome-wise coverage filenames
    depth <- vector("integer", length(coverageFiles))
    cov_files <- list(list())
    
    for (i in seq_along(coverageFiles)) {
        cov_depth <- readRDS(coverageFiles[i])
        depth[i] <- cov_depth$depth
        names(depth)[i] <- names(cov_depth$depth)
        
        tag <- names(cov_depth)[1]
        chr <- names(cov_depth[[1]])
        
        for (seqname in chr){
            cov_files[[seqname]][[tag]] <- cov_depth[[1]][[seqname]]
        }
    }
    
    ## rearrange sample-wise coverage Rle objects into chromosome-wise coverage
    ## Rle objects
    filenames <- list()
    for (seqname in names(cov_files)){
        chr_cov <- list()
        for (tag in names(cov_files[[seqname]])){
            chr_cov[[tag]] <- readRDS(cov_files[[seqname]][[tag]])
        }
        filename <- file.path(od, paste(seqname, "RleCov.RDS"))
        
        ## chr_cov: a list with each of its element containing a Rle of a given
        ## chromosome of each sample
        saveRDS(chr_cov, file = filename)
        filenames[[seqname]] <- filename
    }
    cov_depth <- list(cov = filenames, depth = depth)
    
    ## return a list of filenames pointing to a list of chromosome-wise
    ## coverage of each sample and total read depth of each sample
    cov_depth
}