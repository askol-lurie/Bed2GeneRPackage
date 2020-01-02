library(dplyr)
library(data.table)
library(GenomicRanges)
library(optparse)

getGenes <- function(file){

    ## assumine first column is gene names
    genes <- read.table(file = file, as.is=T, header=F, sep="\t")
    genes <- unlist(c(genes[,1]))
    
    gene_ind <- grep("Gene|gene|GENE", genes)
    if (length(gene_ind) > 0){
        genes <- genes[-gene_ind]
    }

    return(genes)
}


GetGeneInterval <- function(file, keepXtrans = FALSE, keepNR = FALSE){

    ## KEEPXM: SHOULD POSITIONS INCLUDE TRANSCRIPT THAT START WITH XM (COMPUTATIONAL TRANSCRIPTS)

    ## SETTING START AND END TO BE THE MINIMUM AND MAXIMUM POS OBSERVED ##
    ## FOR A GENE
    d <- fread(file, header=T)
    if (keepXtrans == FALSE){
        d <- d %>% filter(substring(name,1,1) != "X")
    }
    if (keepNR == FALSE){
        nrmv <- sum(substring(d$name, 1, 2) == "NR")
        d <- d %>% filter(substring(name,1,2) != "NR")
        print(paste0("Removing ", nrmv, " transcripts starting with NR (non-coding genes)"))
    }
    
    d <- d %>% select(name2, chrom, txStart, txEnd) %>%
        mutate(strand = "+") %>% dplyr::rename(gene = name2) %>%
        group_by(gene, chrom) %>% mutate(start = min(txStart, na.rm=T),
                                         end = max(txEnd, na.rm=T)) %>%
        filter(row_number() == 1) %>% ungroup() %>% select(-txStart, -txEnd)

    ## 

    rng <- makeGRangesFromDataFrame(d, keep.extra.columns = T, seqnames.field = "chrom",
                                      start.field = "start", end.field = "end",
                                   strand.field = "strand")

    ## write.table(rng, file = outFile, quote=F, row.names=F, col.names=T)
    ## print(paste0("Saving gene info to ",outFile))
    return(as.data.frame(rng))
}


gene2bed <- function(genes, geneLocs, prefix, outDir){

    missFile <- paste0(outDir, "/",prefix,"_missingGenes.txt")
    geneIntFile <- paste0(outDir, "/", prefix, ".bed") 
    ## DETERMINE WHICH GENES ARE NOT IN REFSEQ
    GenesMissing <- genes[genes %in% geneLocs$gene == FALSE]

    print(paste0(length(GenesMissing), " genes in gene list but not in RefSeq: "))
    print(paste0("Print list of missing genes to ", missFile))

    ## LIST OF GENES IN MEDEX THAT ARE NOT IN REFSEQ
    write.table(file = missFile, GenesMissing, quote=F, row.names=F, col.names=F)

    ## GET POSITION OF GENES IN GENE LIST
    geneInts <- geneLocs %>% filter(gene %in% genes) %>%
        select(seqnames, start, end, gene) %>% arrange(seqnames, start) %>%
        dplyr::rename(chr = seqnames)
    
    write.table(file = geneIntFile, geneInts, quote=F, row.names=F, col.names=F, sep="\t")
    print(paste0("Wrote intervals for ", nrow(geneInts), " gene to ",geneIntFile))

}
