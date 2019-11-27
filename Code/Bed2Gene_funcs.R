
library(optparse)
library(GenomicRanges)
library(dplyr)

getGenes <- function(geneFiles){

    genes <- c()
    for (file in geneFiles){
        
        tmp <- read.table(file = file, as.is=T, header=FALSE)
        genes <- c(genes, c(unlist(tmp)))
        
    }
    return(genes)
}
        
bed2gene <- function(file, genes, geneLocs, outDir){

    ## FIRST SEE IF THERE ARE ANY MISSING GENES IN GENELOCS ##
    miss <- genes[genes %in% geneLocs$gene == FALSE]
    if (length(miss) > 0){

        fileMiss <- paste0(outDir, "/MissGenes.txt")
        print(paste0("The following genes were not found in the reference gene position list. List of unfound genes printed to ", fileMiss))
        write.table(file = fileMiss, miss, quote=F, row.names = F, col.names=F)
    }
        
    ## SEARCH ONLY FOR THE GENES IN GENES ##
    geneLocs <- geneLocs[which(geneLocs$gene %in% genes)]

    ## GET BED AND CONVERT TO GRANGE OBJECT
    bed <- read.table(file = file, as.is=T, header=F)
    names(bed)[1:3] <- c("chr","start","stop")
    bed$strand = "+"

    bed <- makeGRangesFromDataFrame(bed, keep.extra.columns=FALSE, seqnames.field = "chr", 
                                    start.field = "start", end.field = "stop",
                                    strand.field = "strand")

    ## FIND THE GENE(S) THAT AN INTERVAL OVERLAPS WITH
    ol <- as.data.frame(findOverlaps(bed,geneLocs, type="any", ignore.strand=TRUE ))
    ol <- ol %>% group_by(queryHits) %>%
        mutate(gene = paste(geneLocs$gene[subjectHits],collapse=";")) %>%
        as.data.frame()
    ## ASSIGN GENE(S) TO INTERVALS
    bed$gene = ""
    bed$gene[ol$queryHits] <- ol$gene

    outFile <- basename(file)
    outFile <- gsub("\\..+","_wGenes.bed", outFile)
    outFile <- paste0(outDir,"/", outFile)

    ## DETERMINE IF ANY GENES IN GENES DO NOT HAVE BED INTERVALS ##
    missedGenes <- genes[genes %in% unique(bed$gene) == F]

    if (length(missedGenes) > 0){

        missedGenesTbl <- as.data.frame(geneLocs[which(geneLocs$gene %in% missedGenes)])
        
        missedFile <- basename(file)
        missedFile <- gsub("\\..+","_GenesWoIntervals.txt", missedFile)
        missedFile <- paste0(outDir, "/", missedFile)

        write.table(file = missedFile, missedGenesTbl, quote=F, row.names=F, col.names=T)
        
        print(paste0(length(missedGenes), " gene did not have intervals overlapping them!"))
        print(paste0("Missed genes are written to ",missedFile))
    }
                    
    write.table(file = outFile, bed, quote=F, row.names=F, col.names=F)
    print(paste0("Wrote bed file with gene names to ",outFile))
    return(outFile)
}

makeFastGenePos <- function(file, outFile, keepXtrans = FALSE){

    ## KEEPXM: SHOULD POSITIONS INCLUDE TRANSCRIPT THAT START WITH XM (COMPUTATIONAL TRANSCRIPTS)
    
    ## SETTING START AND END TO BE THE MINIMUM AND MAXIMUM POS OBSERVED ##
    ## FOR A GENE
    d <- fread(file, header=T)
    if (keepXtrans == FALSE){
        d <- d %>% filter(substring(name,1,1) != "X")
    }
    d <- d %>% select(name2, chrom, txStart, txEnd) %>%
        mutate(strand = "+") %>% rename(gene = name2) %>%
        group_by(gene, chrom) %>% mutate(start = min(txStart, na.rm=T),
                                         end = max(txEnd, na.rm=T)) %>%
        filter(row_number() == 1) %>% ungroup() %>% select(-txStart, -txEnd)
    
    ## 
    
    rng <- makeGRangesFromDataFrame(d, keep.extra.columns = T, seqnames.field = "chrom",
                                      start.field = "start", end.field = "end",
                                   strand.field = "strand")

    saveRDS(rng, file = outFile)
    print(paste0("Saving gene info to ",outFile))
}


if (0){

    ## This is how the gene position file was made ##
    
    ## RESOURCE FILES ##
    UCSCGene37fast <- paste0(ResourceDir,"GeneStartEnd37.rds")
    UCSCGeneFile37 <- paste0(ResourceDir,"Genes_GenesPredictions_UCSCRefSeq_GRCh37.gz")

    Gene37 <- makeFastGenePos(UCSCGeneFile37, UCSCGene37fast, keepXtrans = FALSE)
}
