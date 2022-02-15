########################################################
# On Load / Attach
########################################################

InPASDefaults <- list(
  InPAS.logging = TRUE,
  InPAS.genome = NULL,
  InPAS.TxDb = NULL,
  future.globals.maxSize = +Inf,
  InPAS.EnsDb = NULL,
  InPAS.outdir = NULL,
  InPAS.sqliteDb = NULL,
  InPAS.lockname = NULL,
  InPAS.chr2exclude = c(
    "chrM", "MT",
    "chrPltd", "Pltd"
  )
)

#' A function called upon a package is attached to the search path
#'
#' @param libname library name
#' @param pkgname package name
#'
#' @keywords internal
#' @importFrom utils packageVersion
#'
.onAttach <- function(libname, pkgname) {
  # if (!interactive()) return()
  v <- packageVersion("InPAS")
  packageStartupMessage(
    "InPAS : Version ", v,
    "\nFor more information see our website : ",
    "https://bioconductor.org/packages/release/bioc/vignettes/InPAS/inst/doc/InPAS.html\n",
    "If you encounter a bug, please report : https://github.com/jianhong/InPAS/InPAS/issues"
  )
  op <- options()
  toset <- !(names(InPASDefaults) %in% names(op))
  if (any(toset)) {
    options(InPASDefaults[toset])
  }
  invisible(0)
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
#' @param chr2exclude A character vector, NA or NULL, specifying chromosomes or
#'   scaffolds to be excluded for InPAS analysis. `chrM` and alternative scaffolds
#'   representing different haplotypes should be excluded.
#' @import BSgenome
#' @export
addChr2Exclude <- function(chr2exclude = c(
                             "chrM", "MT",
                             "Pltd", "chrPltd"
                           )) {
  if (!is.null(chr2exclude) && !is.character(chr2exclude)) {
    stop("chr2exclude must be NULL, or a character vector!")
  }
  options(InPAS.chr2exclude = chr2exclude)
}

#' Get a globally-applied requirement for filtering scaffolds.
#'
#' This function will get the default requirement of filtering scaffolds.
#'
#' @export
getChr2Exclude <- function() {
  chr2exclude <- options()$InPAS.chr2exclude
  if (is.null(chr2exclude) || is.character(chr2exclude)) {
    return(chr2exclude)
  } else {
    stop(
      "chr to exclude is not set correctly.\n",
      "Please call addChr2Exclude() first"
    )
  }
}

################################################################################
# Add BSgenome as global default
################################################################################

#' Add a globally defined genome to all InPAS functions.
#'
#' This function will set the genome across all InPAS functions.
#'
#' @param genome A BSgenome object indicating the default genome to be
#' used for all InPAS functions. This value is stored as a global environment
#' variable. This can be overwritten on a per-function basis using the given
#' function's `genome` parameter.
#' @import BSgenome
#' @export
addInPASGenome <- function(genome = NULL) {
  if (is(genome, "BSgenome")) {
    message("Setting default genome to ", metadata(genome)$genome, ".")
    options(InPAS.genome = genome)
  } else {
    stop("genome must be a BSgenome object!")
  }
  invisible(0)
}

#' Get the globally defined genome
#'
#' This function will retrieve the genome that is currently in use by InPAS.
#' @export
getInPASGenome <- function() {
  .InPASGenome <- options()$InPAS.genome
  if (!is(.InPASGenome, "BSgenome")) {
    stop(
      "genome is not set. Please call the function\n",
      "addInPASGenome() first"
    )
  }
  .InPASGenome
}

#' Add a globally defined TxDb for InPAS functions.
#'
#' @param TxDb An object of [GenomicFeatures::TxDb-class]
#' @import GenomicFeatures
#' @export
#'
#' @examples
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' addInPASTxDb(TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
addInPASTxDb <- function(TxDb = NULL) {
  if (!is(TxDb, "TxDb")) {
    stop("TxDb must be an object of TxDb!")
  }
  message("Setting default TxDb to ", deparse(substitute(TxDb)), ".")
  options(InPAS.TxDb = TxDb)
  invisible(0)
}

#' Get the globally defined TxDb.
#'
#' @return An object of [GenomicFeatures::TxDb-class]
#' @export
#'
#' @examples
#' library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#' addInPASTxDb(TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene)
#' getInPASTxDb()
getInPASTxDb <- function() {
  TxDb <- options()$InPAS.TxDb
  if (!is(TxDb, "TxDb")) {
    stop(
      "TxDb is not set. Please call the function\n",
      "addInPASTxDb() first"
    )
  }
  TxDb
}

#' Add a globally defined EnsDb to some InPAS functions.
#'
#' @param EnsDb An object of [ensembldb::EnsDb-class]
#'
#' @export
addInPASEnsDb <- function(EnsDb = NULL) {
  if (!is(EnsDb, "EnsDb")) {
    stop("EnsDb must be an object of EnsDb!")
  }
  message("Setting default EnsDb to ", deparse(substitute(EnsDb)), ".")
  options(InPAS.EnsDb = EnsDb)
  invisible(0)
}

#' Get the globally defined EnsDb.
#'
#' @return An object of [ensembldb::EnsDb-class]
#' @export
getInPASEnsDb <- function() {
  EnsDb <- options()$InPAS.EnsDb
  if (!is(EnsDb, "EnsDb")) {
    stop("EnsDb is not set. Please call the function addInPASEnsDb() first")
  }
  EnsDb
}

#' Add a filename for locking a SQLite database
#'
#' @param filename A character(1) vector, specifyong a path to a file for locking.
#' @export
addLockName <- function(filename = NULL) {
  if (is.null(filename) || file.exists(filename)) {
    filename <- tempfile()
    message(
      "The file for locking should not pre-exist!\n",
      "A tempfile ", filename, " is used for locking."
    )
  }
  file.create(filename)
  options(InPAS.lockname = filename)
  invisible(0)
}

#' Get the path to a file for locking the SQLite database
#'
#' @return A path to a file for locking
#' @export

getLockName <- function() {
  filename <- options()$InPAS.lockname
  if (is.null(filename) || !file.exists(filename)) {
    stop(
      "file for locking is not set. please call the function\n",
      "addLockName() first"
    )
  }
  filename
}

#' Add a globally defined output directory to some InPAS functions.
#'
#' @param outdir A character(1) vector, a path with write permission for storing
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @export

addInPASOutputDirectory <- function(outdir = NULL) {
  if (!is.character(outdir) || length(outdir) != 1) {
    stop("An output directory for InPAS analysis is required")
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  outdir <- normalizePath(outdir)
  options(InPAS.outdir = outdir)
}

#' Get the path to a output directory for InPAS analysis
#'
#' @return a normalized path to a output directory for InPAS analysis
#' @export
#'
getInPASOutputDirectory <- function() {
  outdir <- options()$InPAS.outdir
  if (!is.character(outdir) || length(outdir) != 1) {
    stop(
      "outdir is not defined, please call the function\n",
      "addInPASOutputDirectory() first"
    )
  }
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  outdir <- normalizePath(outdir)
  outdir
}

#' Set up global variables for an InPAS analysis
#'
#' @param genome An object [BSgenome::BSgenome-class]. To make things easy, we
#'   suggest users creating a [BSgenome::BSgenome-class] instance from the
#'   reference genome used for read alignment. For details, see the
#'   documentation of [BSgenome::forgeBSgenomeDataPkg()].
#' @param TxDb An object of [GenomicFeatures::TxDb-class]
#' @param EnsDb An object of [ensembldb::EnsDb-class]
#' @param outdir A character(1) vector, a path with write permission for storing
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @param chr2exclude A character vector, NA or NULL, specifying chromosomes or
#'   scaffolds to be excluded for InPAS analysis. `chrM` and alternative scaffolds
#'   representing different haplotypes should be excluded.
#' @param lockfile A character(1) vector, specifying a file name used for
#'   parallel writing to a SQLite database
#'
#' @export
set_globals <- function(genome = NULL,
                        TxDb = NULL,
                        EnsDb = NULL,
                        outdir = NULL,
                        chr2exclude = c(
                          "chrM", "MT",
                          "Pltd", "chrPltd"
                        ),
                        lockfile = tempfile(tmpdir = getInPASOutputDirectory())) {
  addInPASGenome(genome)
  addInPASEnsDb(EnsDb)
  addInPASTxDb(TxDb)
  addInPASOutputDirectory(outdir)
  addLockName(lockfile)
  addChr2Exclude(chr2exclude)
  invisible(0)
}
