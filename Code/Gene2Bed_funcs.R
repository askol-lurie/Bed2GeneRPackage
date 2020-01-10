suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(data.table) )
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(optparse) )

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


GetMergedGeneIntervals <- function(files, keepXtrans = FALSE, keepNR = FALSE){

    d <- c()
    ds <- list()
    for (i in 1:length(files)){        
        ds[[i]] <- GetGeneInterval(files[i], keepXtrans, keepNR)
    }

    for (i in 1:length(ds)){
        
        if (i == 1){
            d <- ds[[i]]
        }else{

            dtmp <- ds[[i]]
            dtmp <- dtmp %>% filter(dtmp$gene %in% d$gene == FALSE)
            d <- rbind(d, dtmp)
        }
    }
               
    return(d)
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
         mutate(strand = "+") %>% dplyr::rename(gene = name2)
    
    ## ADJUST GENES WITH NON OVERLAPPING TRANSCRIPTS. IF A TRANSCRIPT OF A GENE OVERLAPS
    ## WITH NO OTHER TRANSCRIPT OF THAT GENE, THEN PROVIDE IT A NEW NAME:
    ## GENE_X, WHERE X IS AN INTEGER

    d <- adjustGenes(d)
    
    d <- d %>% group_by(gene, chrom) %>% mutate(start = min(txStart, na.rm=T),
                                         end = max(txEnd, na.rm=T)) %>%
        filter(row_number() == 1) %>% ungroup() %>% select(-txStart, -txEnd)

    ## write.table(rng, file = outFile, quote=F, row.names=F, col.names=T)
    ## print(paste0("Saving gene info to ",outFile))
    return(as.data.frame(d))
}

AddTrickyGenes <- function(locs, file, build){

    tg <- read.table(file = file, as.is=T, header=T, sep="\t")
    ind <- grep( paste0("gene|", build), names(tg))
    tg <- tg[,ind]
    names(tg) <- gsub(paste0(build,"_"), "", names(tg))
    tg <- tg %>% mutate(strand = "+", geneExt = gene) %>% select(chrom, strand, gene, geneExt, start, end)
    ind <- which(tg$gene %in% locs$gene == FALSE)
    if (length(ind) > 0){
        locs <- rbind(locs, tg[ind,])
    }
    print(paste0(length(ind), " tricky genes added to list of gene locations."), quote=FALSE)
    return(locs)
}
    
gene2bed <- function(genes, geneLocs, prefix, outDir){

    ## REMOVE ALTERNATIVE LOCI FROM GENELOCS 
    ind <- grep("_", geneLocs$chrom)
    geneLocs <- geneLocs[-ind,]
    
    missFile <- paste0(outDir, "/",prefix,"_missingGenes.txt")
    geneIntFile <- paste0(outDir, "/", prefix, ".bed") 
    ## DETERMINE WHICH GENES ARE NOT IN REFSEQ
    GenesMissing <- genes[genes %in% geneLocs$gene == FALSE]

    print(paste0(length(GenesMissing), " genes in gene list but not in RefSeq."))
    print(paste0("Printing list of missing genes to ", missFile))

    ## LIST OF GENES IN MEDEX THAT ARE NOT IN REFSEQ
    write.table(file = missFile, GenesMissing, quote=F, row.names=F, col.names=F)

    ## GET POSITION OF GENES IN GENE LIST
    geneInts <- geneLocs %>% filter(gene %in% genes) %>%
        mutate(gene = geneExt) %>% 
        select(chrom, start, end, gene) %>% arrange(chrom, start) %>%
        dplyr::rename(chr = chrom)
    
    write.table(file = geneIntFile, geneInts, quote=F, row.names=F, col.names=F, sep="\t")
    print(paste0("Wrote intervals for ", nrow(geneInts), " gene to ",geneIntFile))

}

adjustGenes <- function(data){

    data$row <- 1:nrow(data)
    
    rng <- makeGRangesFromDataFrame(data, keep.extra.columns = TRUE, seqnames.field = "chrom",
                                      start.field = "start", end.field = "end",
                                    strand.field = "strand")

    red <- reduce(rng)
    ol <- findOverlaps(rng, red, type = "any")
    ol <- as.data.frame(ol)

    data <- left_join(data, ol, by=c("row" = "queryHits"))

    data <- data %>% group_by(gene) %>%
        mutate(gap = 1 * (n_distinct(subjectHits) > 1) ,
               int = cumsum(duplicated(subjectHits)==FALSE) ) %>%
        ungroup() %>% group_by(gene, subjectHits) %>%
        mutate(geneExt = ifelse (gap == 0, gene, paste(gene,int,sep="_"))) %>% ungroup()

    data <- data  %>% select(chrom, txStart, txEnd, strand, gene, geneExt)

    return(data)
}
