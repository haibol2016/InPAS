pkgs_required <- c("BSgenome.Mmusculus.UCSC.mm10",
                   "TxDb.Mmusculus.UCSC.mm10.knownGene",
                   "EnsDb.Mmusculus.v79",
                   "rtracklayer",
                   "GenomicRanges",
                   "limma")
for (pkg in pkgs_required)
{
    require(pkg, character.only = TRUE) || stop(pkg, " can't be loaded!")
}

BiocGenerics:::testPackage("InPAS")
