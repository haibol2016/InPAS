#' A package for identifying novel Alternative PolyAdenylation
#' Sites (PAS) based on RNA-seq data
#'
#' The InPAS package provides three categories of important functions:
#' parse_TxDb, extract_UTR3Anno, assemble_allCov, get_ssRleCov, 
#' get_UTR3eSet, test_dPDUI, run_singleSampleAnalysis,
#' run_singleGroupAnalysis, run_limmaAnalysis, filter_testOut,
#' get_usage4plot, setup_GSEA, run_coverageQC
#'
#' @section functions for retrieving 3' UTR annotation: parse_TxDb,
#'   extract_UTR3Anno
#' @section functions for processing read coverage data: assemble_allCov,
#'   get_ssRleCov, run_coverageQC
#' @section  functions for alternative polyadenylation site analysis:
#'    test_dPDUI, run_singleSampleAnalysis, run_singleGroupAnalysis,
#'   run_limmaAnalysis, filter_testOut, get_usage4plot
#'
#' @docType package
#' @name InPAS
globalVariables(c(
  ".", "Predicted_Proximal_APA", "X2",
  "X3", "dup.group", "end.utr3.last",
  "exon_name", "exon_rank", "feature",
  "gene", "start.utr3.last", "transcript",
  "truncated", "tx_name", "Coverage",
  "value"))
