test_get_ssRleCov <- function() {
  ## generate a bedgraph file
  start <- seq(1, 9, by = 2)
  end <- seq(2, 10, by = 2)
  testGR <- GRanges("chr1", IRanges(start, end),
                    strand = "+", score = 1:5)
  filename <- tempfile()
  export(testGR, filename, format = "BEDGraph")
  
  tag <- c("test")
  metadata <- data.frame(tag = tag, 
                         condition = "test",
                         bedgraph_file = filename)
  outdir = tempdir()
  write.table(metadata, file =file.path(outdir, "metadata.txt"), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  sqlite_db <- setup_sqlitedb(metadata = file.path(outdir, "metadata.txt"),
                              outdir)
  
  cov <- get_ssRleCov(bedgraph = filename,
                      tag = "test",
                      genome = BSgenome.Mmusculus.UCSC.mm10,
                      sqlite_db = sqlite_db,
                      outdir = outdir,
                      removeScaffolds = TRUE,
                      BPPARAM = NULL)
  
  checkEquals(rep(2, 5), runLength(cov$test$chr1)[1:5])
  checkEquals(1:5, runValue(cov$test$chr1)[1:5])
}
