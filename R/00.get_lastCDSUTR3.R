#' Extract the last unspliced region of each transcript
#' 
#' Extract the last unspliced region of each transcript from a TxDb. These 
#' regions could be the last 3'UTR exon for transcripts whose 3' UTRs are 
#' composed of multiple exons or last CDS regions and 3'UTRs for transcripts
#' whose 3'UTRs and last CDS regions are on the same single exon.
#'
#' @param TxDb An object of [GenomicFeatures::TxDb-class]
#' @param genome An object of [BSgenome::BSgenome-class]
#' @param chr2exclude A character vector, NA or NULL, specifying chromosomes or 
#'   scaffolds to be excluded for InPAS analysis. `chrM` and alternative scaffolds
#'   representing different haplotypes should be excluded.
#' @param outdir A character(1) vector, a path with write permission for storing 
#'   InPAS analysis results. If it doesn't exist, it will be created.
#' @param MAX_EXONS_GAP An integer(1) vector, maximal gap sizes between the last
#'   known CP sites to a nearest downstream exon. Default is 10 kb for mammalian
#'   genomes. For other species, user need to adjust this parameter.
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @importFrom plyranges as_granges complement_ranges disjoin_ranges filter 
#' group_by mutate reduce_ranges reduce_ranges_directed remove_names select 
#' set_genome_info shift_downstream summarise
#' @importFrom dplyr as_tibble mutate filter arrange bind_rows group_by
#'   left_join summarise n
#' @return A BED file with 6 columns: chr, chrStart, chrEnd, name, score, and
#'   strand.
#' @export

get_lastCDSUTR3 <- function(TxDb = getInPASTxDb(),
                            genome = getInPASGenome(),
                            chr2exclude = getChr2Exclude(),
                            outdir = getInPASOutputDirectory(),
                            MAX_EXONS_GAP = 10000)
{
        if (!is(TxDb, "TxDb")) {
                stop("TxDb must be an object of TxDb!")
        }
        if (!is(genome, "BSgenome"))
        {
                stop("genome must be a BSgenome object!")
        }
        if (!is.null(chr2exclude) && !is.character(chr2exclude))
        {
                stop("chr2Exclude must be NULL or a character vector!")
        }
        if (!is.character(outdir) || length(outdir) != 1) {
                stop(
                        "A path with write permission for storing \n",
                        "the parsed transcriptome is required!"
                )
        }
        if (!dir.exists(outdir)) {
                dir.create(outdir, recursive = TRUE)
        }
        outdir <- normalizePath(outdir)
        seqnames <- trim_seqnames(genome = genome,
                                  chr2exclude = chr2exclude)
        seqlevels <- seqlevels(TxDb)
        seqlevels <- seqlevels[!seqlevels %in% chr2exclude]
        if (!all(seqlevels %in% seqnames))
        {
                stop("chromosome names in TxDb and BSgenome are not consistent!",
                 "\nPlease make sure they match.")
        }
        seqlevels(TxDb) <- seqlevels
        
        ## get all transcripts
        tx <- unlist(transcriptsBy(TxDb, by = "gene")) %>%
                plyranges::mutate(gene = names(.)) %>%
                plyranges::remove_names() %>%
                plyranges::mutate(tx_id = as.character(tx_id)) %>%
                data.frame() %>%
                as_tibble()
        
        ## get all exons and label the last exons
        exons <-
                unlist(exonsBy(TxDb, by = "tx", use.names = TRUE)) %>%
                plyranges::mutate(tx_name = names(.)) %>%
                plyranges::remove_names() %>%
                plyranges::select(-c("exon_id", "exon_name")) %>%
                data.frame() %>%
                as_tibble() %>%
                dplyr::arrange(seqnames, tx_name, -exon_rank)
        
        num_exons <- exons %>%
                dplyr::group_by(tx_name) %>%
                summarise(num_exon = n()) %>%
                as.data.frame()
        rownames(num_exons) <- num_exons[, 1]
        num_exons <- num_exons[, -1, drop = FALSE]
        num_exons <- num_exons[unique(exons$tx_name), , drop = F]
        
        exons$feature <-
                do.call("c", lapply(num_exons[, 1], function(.x) {
                        c("lastutr3", rep("other", time = .x - 1))
                }))
        exons <-
                exons %>% dplyr::left_join(dplyr::select(tx, tx_id, tx_name, gene),
                                           by = "tx_name") %>%
                dplyr::filter(!is.na(gene)) %>%
                dplyr::mutate(exon = paste(seqnames,
                                           feature,
                                           paste(tx_name, exon_rank, sep = "_"),
                                           sep = ":")) %>%
                dplyr::rename(transcript = tx_name) %>%
                plyranges::as_granges()
        
        ## reduce by strand
        last_exons <- exons[exons$feature == "lastutr3"] %>%
                plyranges::reduce_ranges_directed()
        names(last_exons) <-
                paste0("last_exon", seq.int(length(last_exons)))
        
        last_exons_bed <-
                data.frame(
                        chr = as.character(seqnames(last_exons)),
                        start = start(last_exons) - 1,
                        end = end(last_exons),
                        name = paste0("utr", seq.int(length(last_exons))),
                        score = 0,
                        strand = strand(last_exons)
                )
        write.table(
                last_exons_bed,
                file = "last.exons.bed",
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE
        )
        
        ## get next exon gap
        if (is(genome, "BSgenome"))
        {
                genome.info <- as.data.frame(seqinfo(genome))[, 1:2]
        } else {
                genome.info <- as.data.frame(seqinfo(TxDb))[, 1:2]
        }
        genome.info <-
                genome.info[!rownames(genome.info) %in% chr2exclude, ]
        
        gaps <- exons %>%
                plyranges::reduce_ranges() %>%
                plyranges::set_genome_info(
                        seqnames = rownames(genome.info),
                        seqlengths = genome.info$seqlengths,
                        is_circular = genome.info$isCircular
                ) %>%
                plyranges::complement_ranges()
        last_exons.ext1 <- last_exons %>%
                plyranges::shift_downstream(shift = 1L)
        
        ol <- findOverlaps(last_exons.ext1, gaps,
                           ignore.strand = TRUE)
        ol.utr3.clean <- last_exons[queryHits(ol)]
        
        next.exons.gap <- gaps[subjectHits(ol)] %>%
                plyranges::mutate(strand = strand(ol.utr3.clean))
        mcols(next.exons.gap) <- mcols(ol.utr3.clean)
        names(next.exons.gap) <- names(ol.utr3.clean)
        
        wid <- width(next.exons.gap) > MAX_EXONS_GAP
        width(next.exons.gap)[wid &
                                      as.character(strand(next.exons.gap)) == "+"] <-
                MAX_EXONS_GAP
        start(next.exons.gap)[wid &
                                      as.character(strand(next.exons.gap)) == "-"] <-
                end(next.exons.gap)[wid &
                                            as.character(strand(next.exons.gap)) == "-"] - MAX_EXONS_GAP + 1
        next.exons.gap <- next.exons.gap %>%
                plyranges::mutate(feature = "next.exon.gap")
        next.exons.gap.width <- width(next.exons.gap)
        names(next.exons.gap.width) <- names(next.exons.gap)
        next.exons.gap.width <-
                next.exons.gap.width[names(last_exons)]
        next.exons.gap.width <-
                ifelse(is.na(next.exons.gap.width), 0,
                       next.exons.gap.width)
        
        # adjust end for Granges on + strand
        end(last_exons[as.character(strand(next.exons.gap)) == "+"]) <-
                end(last_exons[as.character(strand(next.exons.gap)) == "+"]) +
                next.exons.gap.width[as.character(strand(next.exons.gap)) == "+"] - 1
        # adjust start for Granges on - strand
        start(last_exons[as.character(strand(next.exons.gap)) == "-"]) <-
                start(last_exons[as.character(strand(next.exons.gap)) == "-"]) -
                next.exons.gap.width[as.character(strand(next.exons.gap)) == "-"] + 1
        # reduce
        last_exons <- last_exons %>% plyranges::reduce_ranges()
        last_exons_extended_bed <-
                data.frame(
                        chr = as.character(seqnames(last_exons)),
                        start = start(last_exons) -
                                1,
                        end = end(last_exons),
                        name = paste0("utr", seq.int(length(last_exons))),
                        score = 0,
                        strand = "*"
                )
        write.table(
                last_exons_extended_bed,
                file = "last.exons.extended.reduced.bed",
                sep = "\t",
                quote = FALSE,
                row.names = FALSE,
                col.names = FALSE
        )
        tx_file <-
                file.path(outdir, "00.lastCDS_UTR3.labeled.TxDb.RDS")
        saveRDS(last_exons_extended_bed, file = tx_file)
        last_exons_extended_bed
}
