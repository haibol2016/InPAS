########################################################
# On Load / Attach
########################################################

InPASDefaults <- list(
    InPAS.threads = 1,
    InPAS.logging = TRUE,
    InPAS.genome = NULL,
    InPAS.lockname = tempfile(tmpdir = "."),
    InPAS.chr2exclude = NULL,
    InPAS.verbose = TRUE)

.onAttach <- function(libname, pkgname){
    #if (!interactive()) return()
    v <- packageVersion("InPAS")
    packageStartupMessage("InPAS : Version ", v, 
    "\nFor more information see our website : ",
    "https://bioconductor.org/packages/release/bioc/vignettes/InPAS/inst/doc/InPAS.html\n",
    "If you encounter a bug please report : https://github.com/jianhong/InPAS/InPAS/issues")
    op <- options()
    toset <- !(names(InPASDefaults) %in% names(op))
    if (any(toset)) {options(InPASDefaults[toset])}
    if (!.isWholenumber(options()[["InPAS.threads"]]) || 
       options()[["InPAS.threads"]] == 1 )
    {
        addInPASThreads()
    }
    invisible()
}

#############################################################################
# filtering out scaffolds or not from the genome
#############################################################################
#' Add a globally-applied requirement for filtering out scaffolds from 
#' all analysis
#' 
#' This function will set the default requirement of filtering out scaffolds
#' from all analysis.
#' 
#' @param chr2exclude A character vector indicating the scaffolds/chromosomes to
#' be excluded from all InPAS analyses.
#' @param genome An object of [BSGenome::BSgenome-class]
#' @import BSgenome
#' @export
addChr2Exclude <- function(chr2exclude = NULL, genome = NULL)
{
    if (is.null(genome) || !is(genome, "BSgenome")) 
    {
        stop("genome must be an object of the BSgenome class!")
    }
    if (!is.null(chr2exclude) && 
        (!is.character(chr2exclude) || 
         !any(chr2exclude %in% seqnames(genome)))) 
    {
        stop("chr2exclude must be NULL or valid seqnames compatible with ",
             "seqnames of the genome!")
    }
    options(InPAS.chr2exclude = chr2exclude)
}

#' Get a globally-applied requirement for filtering scaffolds.
#' 
#' This function will get the default requirement of filtering scaffolds.
#' 
#' @export
getChr2Exclude <- function(){
    .chr2Exclude <- options()[["InPAS.chr2exclude"]]
    if (!is.null(.chr2Exclude)){
        if (is.character(.chr2Exclude)){
            .chr2Exclude
        } else {
            NULL
        }
    } else {
        NULL
    }
}

################################################################################
# Parallel Information
################################################################################

#' Add a globally-applied number of threads to use for parallel computing.
#' 
#' This function will set the number of threads to be used for parallel computing
#' across all InPAS functions.
#' 
#' @param threads The default number of threads to be used for parallel execution
#'  across all InPAS functions.
#' This value is stored as a global environment variable.
#' This can be overwritten on a per-function basis using the given function's 
#' `threads` parameter.
#' @param force A logical(1). If you request more than the total number of CPUs
#' minus 2, InPAS will set `threads` to `(nCPU - 2)`. To bypass this, setting 
#' `force = TRUE` will use the number provided to `threads`.
#' @importFrom parallel detectCores
#' @export
addInPASThreads <- function(threads = floor(parallel::detectCores()/ 2),
                            force = FALSE){
    
    if (tolower(.Platform$OS.type) == "windows"){
        message("Detected windows OS, setting threads to 1.")
        threads <- 1
    } else {
        detectedCores <- parallel::detectCores()
        if (threads >= detectedCores - 1){
            if (force) {
                message("Input threads is equal to or greater than ncores minus 1 (",
                        detectedCores - 1,")\nOverriding since force = TRUE.")
            } else {
                message("Input threads is equal to or greater than ncores minus 1 (",
                        detectedCores - 1,")\nSetting cores to ncores minus 2.",
                        "Set force = TRUE to set above this number!")
                threads <- detectedCores - 2
            }
        }
    }
    message("Setting default number of Parallel threads to ", threads, ".")
    options(InPAS.threads = as.integer(round(threads)))
}

#' Get globally-applied number of threads to use for parallel computing.
#' 
#' This function will get the number of threads to be used for parallel execution across all InPAS functions.
#' 
#' @export
getInPASThreads <- function(){
    .InPASThreads <- options()[["InPAS.threads"]]
    if (!is.null(.InPASThreads)){
        if (!.isWholenumber(.InPASThreads)){
            message("option(InPAS.threads) : ", .InPASThreads, 
                    " is not an integer. \nDid you mistakenly set this to a value ",
                    "without addInPASThreads? Resetting to default!")
            addInPASThreads()
            options()[["InPAS.threads"]]
        } else {
            .InPASThreads
        }
    } else {
        1
    }
}

################################################################################
# Add BSgenome as global default
################################################################################

#' Add a globally defined genome to all InPAS functions.
#' 
#' This function will set the genome across all InPAS functions.
#' 
#' @param genome A string or BSgenome object indicating the default genome to be
#' used for all InPAS functions. Currently supported values include "hg19",
#' "hg38","mm9", and "mm10". This value is stored as a global environment 
#' variable. This can be overwritten on a per-function basis using the given 
#' function's `genome` parameter.
#' For something other than one of the currently supported, see `forge_BSgenome`.
#' @param install  A boolean value indicating whether the `BSgenome` object 
#' associated with the provided `genome` should be automatically installed 
#' if it is not currently installed. This is useful for helping reduce
#' user download requirements.
#' @importFrom BiocManager install
#' @export
addInPASGenome <- function(genome = NULL, install = TRUE){
    supportedGenomes <- c("hg19","hg38","mm9","mm10")
    
    if (tolower(genome) %ni% supportedGenomes){
        message("Genome : ", genome, " is not currently supported by InPAS.")
        message("Currently supported genomes : ", paste0(supportedGenomes, 
                                                         collapse = ","))
        if (is(genome, "BSgenome"))
        {
            message("Setting default genome to ", metadata(genome)$genome, ".")
            options(InPAS.genome = genome)
        } else {
            stop("The genome is neither currently supported nor an object of BSgenome. ",
            "To continue try building a custom BSgenome with forge_BSgenome() ",
            "and install it first!")
        }
    } else {
        #Check if BSgenome exists!
        if (tolower(genome)=="hg19"){
            if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg19", 
                                 quietly = TRUE)){
                if (install){
                    message("BSgenome for hg19 not installed! Now installing by ",
                            "the following:\n\tBiocManager::install",
                            "(\"BSgenome.Hsapiens.UCSC.hg19\")")
                    BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
                } else {
                    stop("BSgenome for hg19 not installed! Please install by ",
                         "setting install = TRUE or by the following:\n\t",
                         "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg19\")")
                }
            }
        } else if (tolower(genome)=="hg38"){
            if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", 
                                 quietly = TRUE)){
                if (install){
                    message("BSgenome for hg38 not installed! Now installing by ",
                            "the following:\n\tBiocManager::install", 
                            "(\"BSgenome.Hsapiens.UCSC.hg38\")")
                    BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
                } else {
                    stop("BSgenome for hg38 not installed! Please install by ",
                         "setting install = TRUE or by the following:\n\t", 
                         "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")")
                }
            }
        } else if (tolower(genome)=="mm9"){
            if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm9", 
                                 quietly = TRUE)){
                if (install){
                    message("BSgenome for mm9 not installed! Now installing by ",
                            "the following:\n\tBiocManager::install", 
                            "(\"BSgenome.Mmusculus.UCSC.mm9\")")
                    BiocManager::install("BSgenome.Mmusculus.UCSC.mm9")
                } else {
                    stop("BSgenome for mm9 not installed! Please install by ",
                         "setting install = TRUE or by the following:\n\t",
                         "BiocManager::install(\"BSgenome.Mmusculus.UCSC.mm9\")")
                }
            }
        } else if (tolower(genome)=="mm10"){
            if (!requireNamespace("BSgenome.Mmusculus.UCSC.mm10",
                                 quietly = TRUE)){
                if (install){
                    message("BSgenome for mm10 not installed! Now installing by ",
                            "the following:\n\tBiocManager::install",
                            "(\"BSgenome.Mmusculus.UCSC.mm10\")")
                    BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
                } else {
                    stop("BSgenome for mm10 not installed! Please install by ",
                         "setting install = TRUE or by the following:\n\t",
                         "BiocManager::install(\"BSgenome.Mmusculus.UCSC.mm10\")")
                }
            }
        }
        genome <- paste0(toupper(substr(genome, 1, 1)), 
                         substr(genome, 2, nchar(genome)))
        message("Setting default genome to ", genome, ".")
        options(InPAS.genome = genome)
    }
    invisible(0)
}

#' Get the globally defined genome
#' 
#' This function will retrieve the genome that is currently in use by InPAS. 
#' @export
getInPASGenome <- function(){
    supportedGenomes <- c("hg19","hg38","mm9","mm10")
    .InPASGenome <- options()[["InPAS.genome"]]
    
    if (!is.null(.InPASGenome)){
        ag <- .InPASGenome
        if (!is.character(ag) && !is(ag, "BSgenome")){
            return(NULL)
        } else if (is.character(ag)) {
             if(tolower(ag) %in% supportedGenomes){
                genome <- paste0(toupper(substr(ag, 1, 1)), 
                                 substr(ag, 2, nchar(ag)))
                return(genome)
            } else {
                stop("option(InPAS.genome) : ", ag, 
                     " is not currently supported by InPAS.", 
                     "\nDid you mistakenly set this to a value without ",
                     "addInPASGenome()?")
            }
        } else {
            return(ag)
        }
    } else {
        return(NULL)
    }
}


#' Add a filename for locking
#'
#' @param filename A character(1) vector, a path to a file for locking.
#' @export
addLockName <- function(filename = NULL)
{
    if (is.null(filename) || file.exists(filename))
    {
        filename <- tempfile()
        message("The file for locking should not pre-exist!\n",
                "A tempfile ", filename, " is used for locking.")
    } else {
        file.create(filename)
    }
    options()["InPAS.lockname"] <- filename
    invisible(0)
}


#' Get the path to a file for locking
#'
#' @return A path to a file for locking
#' @export

getLockName <- function()
{
    options()["InPAS.lockname"]
}

