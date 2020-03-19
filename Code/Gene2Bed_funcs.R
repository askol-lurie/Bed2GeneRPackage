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


GetMergedGeneIntervals <- function(files, keepXtrans = FALSE, keepNR = TRUE, coding=FALSE){

    d <- c()
    ds <- list()
    genes <- c()
    for (i in 1:length(files)){        
        ds[[i]] <- GetGeneInterval(files[i], keepXtrans, keepNR, coding, genes)
        genes <- c(genes, unique(ds[[i]]$gene))
    }

    for (i in 1:length(ds)){
        
        if (i == 1){
            d <- ds[[i]]
        }else{
            if (nrow(d) > 0){ ## make sure that new genes are added by file
                dtmp <- ds[[i]]
                dtmp <- dtmp %>% filter(dtmp$gene %in% d$gene == FALSE)
                d <- rbind(d, dtmp)
            }
        }
    }
               
    return(d)
}

GetGeneInterval <- function(file, keepXtrans = FALSE, keepNR = FALSE,
                            coding=FALSE, genes){

    ## KEEPXM: SHOULD POSITIONS INCLUDE TRANSCRIPT THAT START WITH XM (COMPUTATIONAL TRANSCRIPTS)

    ## coding == false: REGION RETURNED IS TRANSCRIPTION START AND END (INCLUDES UTRS)
    ## coding == true: REGION RETURNED IS CODING START AND END
    
    ## SETTING START AND END TO BE THE MINIMUM AND MAXIMUM POS OBSERVED ##
    ## FOR A GENE
    print(paste0("Processing file ",file), quote=FALSE)
    d <- fread(file, header=T)

    ## remove genes that have already been seen in previous files since
    ## files are processed in order of preference
    ## name2 = gene in d
    d <- d %>% filter(name2 %in% genes == FALSE)

    ## make sure there are additional genes to be processed (e.g. not all genes
    ## removed
    no.new.trans <- nrow(d)
    print(paste0("Number of new transcripts to add: ",no.new.trans), quote=FALSE)
    if (no.new.trans == 0){

        d = c()
    }else{
        
        if (keepXtrans == FALSE){
            d <- d %>% filter(substring(name,1,1) != "X")
        }
        if (keepNR == FALSE){
            nrmv <- sum(substring(d$name, 1, 2) == "NR")
            d <- d %>% filter(substring(name,1,2) != "NR")
            print(paste0("Removing ", nrmv, " transcripts starting with NR (non-coding genes)"))
        }
        
        if (coding == TRUE){
            d <- d %>% select(name2, chrom, cdsStart, cdsEnd, strand) %>%
                dplyr::rename(gene = name2, start = cdsStart, end = cdsEnd)
        }else{
            d <- d %>% select(name2, chrom, txStart, txEnd, strand)  %>%
                dplyr::rename(gene = name2, start = txStart, end = txEnd)
        }
        
        ## ADJUST GENES WITH NON OVERLAPPING TRANSCRIPTS. IF A TRANSCRIPT OF A GENE OVERLAPS
        ## WITH NO OTHER TRANSCRIPT OF THAT GENE, THEN PROVIDE IT A NEW NAME:
        ## GENE_X, WHERE X IS AN INTEGER
        
        d <- adjustGenes(d)
        
        d <- d %>% group_by(geneExt, chrom) %>%
            mutate(startMin = min(start, na.rm=T),
                   endMax = max(end, na.rm=T)) %>%
            filter(row_number() == 1) %>% ungroup() %>% select(-start, -end) %>%
            dplyr::rename(start = startMin, end = endMax)
    }
    return(as.data.frame(d))
}

GetCodingIntervals <- function(files, keepXtrans = FALSE, keepNR = FALSE){

    d <- c()
    ds <- list()
    for (i in 1:length(files)){        
        ds[[i]] <- GetCodingInterval(files[i], keepXtrans, keepNR)
    }

    for (i in 1:length(ds)){
        
        if (i == 1){
            d <- ds[[i]]
        }else{

            dtmp <- ds[[i]]
            ## assuming that files are ordered by importance.
            ## if a gene is in a previous file, don't include additional coding intervals
            ## from current file
            dtmp <- dtmp %>% filter(dtmp$gene %in% d$gene == FALSE)
            d <- rbind(d, dtmp)
        }
    }

    ## merge intervals
    d <- reduce(d)
    
    return(d)
}

GetCodingInterval <- function(file, keepXtrans = FALSE, keepNR = FALSE){

    ## KEEPXM: SHOULD POSITIONS INCLUDE TRANSCRIPT THAT START WITH XM (COMPUTATIONAL TRANSCRIPTS)
    
    ## SETTING START AND END TO BE THE MINIMUM AND MAXIMUM POS OBSERVED ##
    ## FOR A GENE
    d <- fread(file, header=T)
    if (keepXtrans == FALSE){
        nrmv <- sum(substring(d$name,1,1) == "X")
        d <- d %>% filter(substring(name,1,1) != "X")
        print(paste0("Removing ", nrmv, " transcripts starting with X (predicted genes)"))
        
    }
    if (keepNR == FALSE){
        nrmv <- sum(substring(d$name, 1, 2) == "NR")
        d <- d %>% filter(substring(name,1,2) != "NR")
        print(paste0("Removing ", nrmv, " transcripts starting with NR (non-coding genes)"))
    }

    ## change d so that each row is an exon
    de <- codeify(d)
    
    rng <- makeGRangesFromDataFrame(de, keep.extra.columns = T, seqnames.field = "chr",
                                    start.field = "start", end.field = "end",
                                    strand.field = "strand")

    return(rng)
}


codeify <- function(data){

    de <- list()
    for (i in 1:nrow(data)){
        if (i%%5000 == 0){
            print(paste0("Working on transcript ",i, " of ",nrow(data)))
        }
        tran = data$name[i]
        gene = data$name2[i]
        strand = data$strand[i]
        chrom = data$chrom[i]
        cdsStart = data$cdsStart[i]
        cdsEnd = data$cdsEnd[i]
        exonFrames <- strsplit(data$exonFrames[i], split=",")[[1]]
        tmp <- cbind(strsplit(data$exonStarts[i], split=",")[[1]],
                     strsplit(data$exonEnds[i], split=",")[[1]])

        ## remove exons that are UTRs only ##
        rmInd <- which(exonFrames == -1)
        ## if there is no coding region don't do anything
        if (length(rmInd) != length(exonFrames)){
            ## if there are exons that are all UTR, remove them ##
            if (length(rmInd) > 0){
                tmp <- tmp[-rmInd,]
            }
            ## if only one coding region remains
            if (length(tmp) == 2){
                tmp <- matrix(c(cdsStart, cdsEnd), 1, 2)
            }else{
                tmp[1,1] <- cdsStart
                tmp[nrow(tmp),2] <- cdsEnd
            }
                       
            de[[i]] <- cbind(chrom, tmp, gene, strand, tran)
        }
    }
    
    de <- do.call(rbind, de)
    colnames(de) <- c("chr","start","end","gene","strand","tran")
    de <- as.data.frame(de, string.as.factors = FALSE)    
    
    ## keep distinct intervals (means will delete transcript names) ##
    print("Removing duplicate intervals (as a result of mult transcripts per gene)")
    de <- de %>% group_by(chr,start,end,gene) %>% filter(row_number() == 1) %>% ungroup()
    
    return(de)
}
      
    
gene2bed <- function(genes, geneLocs, prefix, outDir, pad=0, rmChrM=TRUE){

    ## REDUCE GENES TO UNIQUE GENES AND REMOVE NAS OR "" ##
    genes <- unique(na.omit(genes))
    ind <- which(genes == "")
    if (length(ind) > 0){ genes <- genes[-ind] }
    
    ## REMOVE ALTERNATIVE LOCI FROM GENELOCS
    ind <- grep("_", geneLocs$chrom)
    keepMito <- grep("NC_", geneLocs$chrom)
    ind <- ind[ind %in% keepMito == FALSE]

    ## REMOVE GENES WITH CHRM (USING NC_012920) INSTEAD
    if (rmChrM == TRUE){
        ind <- unique(c(ind, which(geneLocs$chrom == "chrM")))
    }
    
    print(paste0("Removing ", length(ind), " genes on alternative and chrM (NC_012920 used instead)"))
    
    genesRm <- c()
    if (length(ind) > 0){
        genesRm <- unique(geneLocs$gene[ind])
        geneLocs <- geneLocs[-ind,]
        
        ## IF ANY OF THE REMAINING GENES HAD ALT LOC VERSIONS REMOVED THEY MAY HAVE AN
        ## UNDERSCORE THAT CAN BE REMOVED
        geneLocs <- geneLocs %>% group_by(gene) %>%
            mutate(nints = n_distinct(geneExt),
                   geneExt = ifelse(nints == 1, gene, geneExt)) %>%
            ungroup() %>% as.data.frame()
    }

    ## REPORT IF ANY GENES ARE REMOVED THAT ARE ONLY ON ALTERNATIVE LOCI
    ind <- which( (genes %in% genesRm) & (genes %in% geneLocs$gene == FALSE) )

    if (length(ind) > 0){
        altLocFile <- paste0(outDir, "/", prefix,"_altLocGenes.txt")
        print(paste0(length(ind), " genes are on alternative loci only."))
        write.table(file = altLocFile, genes[ind], quote=F, row.names=F, col.names=F)
    }

    ## DETERMINE WHICH GENES ARE NOT IN REFSEQ
    missFile <- paste0(outDir, "/",prefix,"_missingGenes.txt")
    geneIntFile <- paste0(outDir, "/", prefix, ".bed") 

    GenesMissing <- genes[genes %in% geneLocs$gene == FALSE]

    print(paste0(length(GenesMissing), " genes in gene list but not in RefSeq."))
    print(paste0("Printing list of missing genes to ", missFile))

    ## LIST OF GENES IN GENE LIST THAT ARE NOT IN REFSEQ
    write.table(file = missFile, GenesMissing, quote=F, row.names=F, col.names=F)

    ## GET POSITION OF GENES IN GENE LIST
    geneInts <- geneLocs %>% filter(gene %in% genes) %>%
        mutate(gene = geneExt) %>% 
        select(chrom, start, end, gene) %>% arrange(chrom, as.numeric(start)) %>%
        dplyr::rename(chr = chrom)

    ## ADD PADDING IF REQUESTED ##
    if (pad != 0){
        print(paste0("Adding ",pad, " bp padding to gene start and end positions"))
        geneInts <- geneInts %>% mutate(start = start - pad, end = end + pad)
    }
    
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

    data <- data  %>% select(chrom, start, end, strand, gene, geneExt)

    return(data)
}


makeGeneLocFile <- function(files, mitoFile = "", ResourceDir, Prefix,
                            build, keepXtrans = FALSE, keepNR = TRUE, coding = FALSE){

    outFile <- paste0(ResourceDir,Prefix,"_",build,".rds")
    geneLocs <- GetMergedGeneIntervals(files, keepXtrans, keepNR, coding)

    if (mitoFile != ""){
        
        print(paste0("Adding mitochondrial genes from ",mitoFile))
        mito <- read.table(file = mitoFile, as.is=T, header=FALSE)
        names(mito) <- c("chrom", "start","end","gene")
        mito$strand <- "+"
        mito$geneExt <- mito$gene
        mito <- mito %>% select(chrom, strand, gene, geneExt, start, end)

        geneLocs <- rbind(geneLocs, mito)
      }
    
    saveRDS(file = outFile, geneLocs)

    print(paste("Wrote gene location file to ",outFile))
}


makeCodingLocFile <- function(files, ResourceDir, Prefix, build, keepXtrans = FALSE, keepNR = TRUE){

    outFile <- paste0(ResourceDir,Prefix,"_",build,".rds")
    geneLocs <- GetCodingIntervals(files, keepXtrans, keepNR)
    saveRDS(file = outFile, geneLocs)

    print(paste("Wrote gene location file to ",outFile))
}

getCodingLocs <- function(geneCodingFile){

    c <- readRDS(geneCodingFile)

    c <- c %>% select(chr, strand, gene, startCd, endCd) %>%
        dplyr::rename(chrom = chr, start = startCd, end = endCd) %>%
        filter(is.na(start) == FALSE) %>%
        group_by(gene) %>% mutate(row = row_number()) %>% ungroup() %>%
        mutate(geneExt = paste0(gene, "_", row)) %>%
        select(chrom, strand, gene, geneExt, start, end)

    return(c)
}
