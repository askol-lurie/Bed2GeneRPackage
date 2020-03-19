
suppressPackageStartupMessages( library(optparse) ) 
suppressPackageStartupMessages( library(GenomicRanges) )
suppressPackageStartupMessages( library(dplyr) )
suppressPackageStartupMessages( library(data.table) )



getGenes <- function(geneFiles){

    genes <- c()
    for (file in geneFiles){
        
        tmp <- read.table(file = file, as.is=T, header=FALSE)
        genes <- c(genes, c(unlist(tmp)))
        
    }
    return(genes)
}
        
bed2gene <- function(file, genes = c(), codingLocs, geneLocs,
                     prefix = "", outDir, rmChrM=TRUE){

    ## SET UP OUTPUT FILE PREFIX ##
    if (prefix == ""){
        outFile <- basename(file)        
        outFile <- gsub("\\..+","", outFile)
    }else{
        outFile <- prefix
    }
    outFile <- paste0(outFile , "_")

    ## REMOVE ALTERNATIVE LOCI FROM GENELOCS
    ind <- grep("_", codingLocs$chr)
    keepMito <- grep("NC_", codingLocs$chr)
    ind <- ind[ind %in% keepMito == FALSE]

    ## REMOVE GENES WITH CHRM (USING NC_012920) INSTEAD
    if (rmChrM == TRUE){
        ind <- unique(c(ind, which(codingLocs$chr == "chrM")))
    }
    
    print(paste0("Removing ", length(ind), " genes on alternative and chrM (NC_012920 used instead) from gene location file."))
    
    ## IF GENE LIST PROVIDED THEN
    ## 1. TRIM GENELOCS TO ONLY INCLUDE THOSE GENES,
    ## 2. DETERMINE IF ANY GENES IN GENE LIST AREN'T IN GENELOCS,
    ## 3. CREATE BED FILE CONTAIN THE EXONIC REGIONS OF GENES IN GENELIST
    if (length(genes) > 0){
        
        ## FIRST SEE IF THERE ARE ANY MISSING GENES IN GENELOCS ##
        miss <- genes[genes %in% codingLocs$gene == FALSE]
        if (length(miss) > 0){
            
            fileMiss <- paste0(outDir, "/MissGenes.txt")
            print(paste0("The following genes were not found in the reference gene position list :"))
            print(miss)
            print(paste0("List of unfound genes printed to ", fileMiss))
            write.table(file = fileMiss, miss, quote=F, row.names = F, col.names=F)
        }
        
        ## SEARCH ONLY FOR THE GENES IN GENES ##
        codingLocs <- codingLocs[which(codingLocs$gene %in% genes), ]
        geneLocs <- geneLocs[which(geneLocs$gene %in% genes), ]
        
        ## WRITE OUT A BED FILE BASED ON THE EXONIC INTERVALS
        if (0){  ## not sure why this is needed. Delete?
            codingLocBed <- codingLocs %>%
                select(chr, start, end, strand, gene, tran) %>%
                dplyr::rename(stop = end)
            
            codingLocBed <- codingLocBed %>% group_by(chr, start, stop, gene) %>%
                filter(row_number() == 1) %>%
                select(chr, start, stop, strand, gene)
            FileOut = paste0(outDir,"/", outFile,"ExonIntervals.bed")
            write.table(file = FileOut, codingLocBed, quote=F, row.names=F,
                        col.names=TRUE, sep="\t")
        }
    }
    
    ## GET BED AND CONVERT TO GRANGE OBJECT
    bed <- read.table(file = file, as.is=T, header=F)
    names(bed)[1:3] <- c("chr","start","stop")
    bed$strand = "+"

    ## Remove any extra columns that were in the bed file ##
    bed <- bed %>% select(chr, start, stop, strand)
    
    bed <- makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE,
                                    seqnames.field = "chr", 
                                    start.field = "start", end.field = "stop",
                                    strand.field = "strand")
    codingLocs <- makeGRangesFromDataFrame(codingLocs, keep.extra.columns=TRUE,
                                           seqnames.field = "chr", 
                                           start.field = "start",
                                           end.field = "end",
                                           strand.field = "strand")
    geneLocs <- makeGRangesFromDataFrame(geneLocs, keep.extra.columns=TRUE,
                                           seqnames.field = "chrom", 
                                           start.field = "start",
                                           end.field = "end",
                                         strand.field = "strand")
    ##
    ## FIND THE CODING REGION(S) THAT AN INTERVAL OVERLAPS WITH
    ##
    ol <- as.data.frame(findOverlaps(bed,codingLocs, type="any",
                                     ignore.strand=TRUE ))
    o.bp <- as.data.frame( pintersect(bed[ol[,1],], codingLocs[ol[,2],]) )
    ol <- ol %>% mutate(ol.bp = o.bp$width)
    ol <- ol %>% mutate(gene = codingLocs$gene[subjectHits])

    ## remove multiple lines for each unique inverval-gene entry
    ol <- ol %>% group_by(queryHits, gene) %>%
        mutate(ol.bp = max(ol.bp)) %>% filter(row_number() == 1)

    ## if multiple genes map to the same interval then combined gene names and
    ## remove multiple occurance of interval
    ol <- ol %>% group_by(queryHits) %>%
        mutate(gene = ifelse(n() > 1, paste(gene,collapse=";"),
                             gene),
               ol.bp = as.character(ol.bp),
               ol.bp = ifelse(n() > 1, paste(ol.bp, collapse=";"),
                              ol.bp))%>%
        filter(row_number() == 1)

    ##
    ## FIND THE GENE REGION(S) THAT AN INTERVAL OVERLAPS WITH
    ##
    ol.gene <- as.data.frame(findOverlaps(bed, geneLocs, type="any",
                                     ignore.strand=TRUE ))
    o.bp <- as.data.frame( pintersect(bed[ol.gene[,1],],
                                      geneLocs[ol.gene[,2],]) )
    ol.gene <- ol.gene %>% mutate(ol.bp = o.bp$width)
    ol.gene <- ol.gene %>% mutate(gene = geneLocs$gene[subjectHits])

    ## remove multiple lines for each unique inverval-gene entry
    ol.gene <- ol.gene %>% group_by(queryHits, gene) %>%
        mutate(ol.bp = max(ol.bp)) %>% filter(row_number() == 1)

    ## if multiple genes map to the same interval then combined gene names and
    ## remove multiple occurance of interval
    ol.gene <- ol.gene %>% group_by(queryHits) %>%
        mutate( gene = ifelse(n() > 1, paste(gene,collapse=";"), gene) ,
               ol.bp = as.character(ol.bp), 
               ol.bp = ifelse(n() > 1, paste(ol.bp, collapse=";"), ol.bp) )%>%
        filter(row_number() == 1)
    
    ## ASSIGN CODING GENE(S) TO INTERVALS
    bed$gene_coding <- ""
    bed$gene_coding[ol$queryHits] <- ol$gene

    ## ASSIGN GENES TO INTERVALS
    bed$gene_txn <- ""
    bed$gene_txn[ol.gene$queryHits] <- ol.gene$gene
    
    ## CHANGE CHROMOSOME COLUMN NAME BACK TO CHR INSTEAD OF SEQNAMES ##
    bed <- bed %>% as.data.frame() %>% rename(chr = seqnames)
    
    FileOut <- paste0(outDir,"/", outFile, "AllInts_wGenes.bed")
    write.table(file = FileOut, bed, quote=F, row.names=F, col.names=T,
                sep="\t")
    print(paste0("Wrote bed file with gene names to ",FileOut))

    print(paste0("***** number of genes in genes ",length(genes)))
    if (length(genes) > 0){
        ## DETERMINE IF ANY GENES IN GENE LIST DO NOT HAVE BED INTERVALS ##
        missedGenesCoding <- genes[genes %in% unique(bed$gene_coding) == F]
        missedGenesTxn <- genes[genes %in% unique(bed$gene_txn) == F]
        
        if (length(missedGenes) > 0){

            ## REPORT GENES IN GENE LIST WITH NO CODING REGIONS
            ## IN BED FILE
            missedGenesTbl <- as.data.frame(codingLocs[which(codingLocs$gene %in% missedGenesCoding)])
            
            FileOut <- paste0(outDir,"/", outFile, "GenesCodingWoIntervals.txt")
            write.table(file = FileOut, missedGenesTbl, quote=F, row.names=F, col.names=T, sep="\t")
            
            print(paste0(length(missedGenes), " gene coding regions did not have intervals overlapping them!"))
            print(paste0(unique(missedGenesTbl$gene), collapse=", "))
            print(paste0("Missed genes coding regions are written to ",FileOut))

            ## REPORT GENES IN GENE LIST WITH REGIONS
            ## IN BED FILE (ANYWHERE IN GENE)
            missedGenesTbl <- as.data.frame(codingLocs[which(codingLocs$gene %in% missedGenesTxn)])
            
            FileOut <- paste0(outDir,"/", outFile, "GenesTxnWoIntervals.txt")
            write.table(file = FileOut, missedGenesTbl, quote=F, row.names=F, col.names=T, sep="\t")
            
            print(paste0(length(missedGenes), " gene did not have intervals overlapping them!"))
            print(paste0(unique(missedGenesTbl$gene), collapse=", "))
            print(paste0("Missed genes are written to ",FileOut))
        }                
    }

    ## convert start and end sites to interger to ensure non-scientific notation
    bed <- bed %>% mutate(start = as.integer(start),
                          end = as.integer(end) )
    
    ## IF BED FILE CONTAINED GENES (COLUMN NAME gene)
    ## then also output bed intervals with those genes (given name) only ##
    bed <- bed %>% as.data.frame()
    if (any(names(bed) == "gene")){

        bed <- bed %>% filter(gene_txn %in% genes)
        FileOut <- paste0(outDir,"/", outFile, "GivenGenesOnly.bed")
        write.table(file = FileOut, bed, quote=F, row.names=F, col.names=T, sep="\t")
        print("Writing bed intervals containing genes from gene file")
        print(paste0("Wrote ",FileOut))
    }    

    return()
}



makeExonLocFile <- function(files, ResourceDir, mitoFile = "",
                            Prefix = "GeneExons", build, keepXtrans = FALSE,
                            keepNR = TRUE){
    
    ## This is how the gene position file was made ##

    Exons <- c()
    genes <- c() ## keep track of genes so that genes that have already
                 ## been seen don't have to be processed 
    for (i in 1:length(files)){

        file <- files[i]
        print(paste0("Processing file ",file), quote=F)

        source <- basename(file)
        source <- gsub("\\..+","",source)
        
        tmp <- makeFastGenePos(file, keepXtrans = keepXtrans, keepNR = keepNR,
                               genes = genes)
        
        if (i != 1 & is.null(tmp) == FALSE){
            ind <- which(tmp$gene %in% Exons$gene == FALSE)
            ##if (length(ind) > 0){
            ##    tmp <- tmp[ind,]
            tmp$source <- source
            Exons <- rbind(Exons, tmp)
                
            genes <- c(genes, unique(tmp$gene))
        }        
    }
    
    if (mitoFile != ""){
        print(paste0("Adding mitochondrial genes from ",mitoFile))
        mito <- read.table(file = mitoFile, as.is=T, header=FALSE)
        names(mito) <- c("chr", "start","end","gene")
        mito$startCd = mito$start
        mito$endCd <- mito$end
        mito$exon <- 1
        mito$strand <- "+"
        mito$tran <- mito$gene
        mito <- mito %>% select(chr, start, end, startCd, endCd, exon, gene, strand, tran)
        source <- basename(mitoFile)
        source <- gsub("\\..+","",source)
        mito$source <- source
        Exons <- rbind(Exons, mito)
    }
    
    outFile <- paste0(ResourceDir, Prefix,"_",build,".rds")
    saveRDS(file = outFile, Exons)
    
    print(paste0("Created exon resource files ",outFile), quote=FALSE)    
}

makeFastGenePos <- function(file, outFile, keepXtrans = FALSE,
                            keepNR = FALSE, genes = genes){
        
    ## KEEPXM: SHOULD POSITIONS INCLUDE TRANSCRIPT THAT
    ## START WITH XM (COMPUTATIONAL TRANSCRIPTS)
    
    ## SETTING START AND END TO BE THE MINIMUM AND MAXIMUM POS OBSERVED ##
    ## FOR A GENE
    d <- fread(file, header=T)
    ## name2 becomes gene
    ## remove genes that have been seen in previous files, since info in
    ## previous files is prefered
    d <- d %>% filter(d$name2 %in% genes == FALSE)

    no.new.trans <- nrow(d)
    print(paste0("Additional transcripts being added: ",no.new.trans),
          quote=FALSE)    
    if (no.new.trans == 0){
        de <- c()
    }else{
        
        if (keepXtrans == FALSE){
            nrmv <- sum(substring(d$name,1,1) == "X")
            d <- d %>% filter(substring(name,1,1) != "X")
            print(paste0("Removing ", nrmv, " transcripts starting with X (predicted genes)"), quote=FALSE)
            
        }
        if (keepNR == FALSE){
            nrmv <- sum(substring(d$name, 1, 2) == "NR")
            d <- d %>% filter(substring(name,1,2) != "NR")
            print(paste0("Removing ", nrmv, " transcripts starting with NR (non-coding genes)"), quote=F)
        }
        
        ## change d so that each row is an exon
        de <- exonify(d)
    }
    return(de)
}

exonify <- function(data){

    de <- list()
    for (i in 1:nrow(data)){
        if (i%%5000 == 0){
            print(paste0("Working on transcript ",i, " of ",nrow(data)))
        }
        tran <- data$name[i]
        gene <- data$name2[i]
        strand <-  data$strand[i]
        chrom <- data$chrom[i]
        cdStart <- data$cdsStart[i]
        cdEnd <- data$cdsEnd[i]
        exonPos <- cbind(strsplit(data$exonStarts[i], split=",")[[1]],
                     strsplit(data$exonEnds[i], split=",")[[1]])
        frames <- strsplit(data$exonFrames[i], split=",")[[1]]
        exonNo <- 1:data$exonCount[i]
        if (strand == "-"){exonNo <- rev(exonNo)}
        cdPos <- exonPos
        cdInd <- which(frames != -1)
        if (length(cdInd) < data$exonCount[i]){
            cdPos[-cdInd,] <- NA
        }
        if (length(cdInd) > 0){
            cdStartInd <- min(cdInd)
            cdEndInd <- max(cdInd)
            cdPos[cdStartInd,1] <- cdStart
            cdPos[cdEndInd,2] <- cdEnd
        }
            
        de[[i]] <- cbind(chrom, exonPos, cdPos, exonNo, gene, strand, tran)
    }
    
    de <- do.call(rbind, de)

    de <- as.data.frame(de, string.as.factors = FALSE)    
    names(de) <- c("chr","start","end","startCd","endCd", "exon", "gene","strand","tran")

    ## BEFORE REMOVING DISTINCT INTERALS, KEEP TRACK OF THE TRANSCRIPT NAMES IN WHICH THEY ARE
    ## EXONS AND WHICH EXONS THEY ARE
    de <- de %>% mutate(tran = as.character(tran),
                        exon = as.character(exon))
    ##
    print("Collating transcript ids and exons. . .", quote=F)
    deMultTrans <- de %>% group_by(gene, chr,start,end) %>% filter(n_distinct(tran) > 1) %>%
        select(chr, start, end, gene, exon, tran)
    deMultTrans <- deMultTrans %>%
        summarize(tranCt = paste(tran, collapse=";"), exonCt = paste(exon, collapse=";"))
    deMultTrans <- deMultTrans %>% ungroup()
    print("Done collating.", quote=F)
    
    ## Replace concatinates transcripts and exons for
    de$startCd = as.numeric(as.character(de$startCd))
    de$endCd = as.numeric(as.character(de$endCd))
    de <- de %>% group_by(chr,start,end,gene) %>%
        summarize(startCd = ifelse(all(is.na(startCd)) , NA, min(startCd, na.rm=T)),
                  endCd   = ifelse(all(is.na(endCd)),    NA, max(endCd, na.rm=T)),
                  tran = tran[1], strand = strand[1], exon=exon[1]) %>% ungroup()
                                                      
    de <- left_join(de, deMultTrans, by = c("chr", "start", "end", "gene"))
    de <- de %>% mutate(tran = ifelse(!is.na(tranCt), tranCt, tran),
                        exon = ifelse(!is.na(exonCt), exonCt, exon) ) %>%
        select(-tranCt, -exonCt)
    
    ## keep distinct intervals (means will delete transcript names) ##
    print("Removing duplicate intervals (as a result of mult transcripts per gene)")
    de <- de %>% mutate(chr = as.character(chr),
                        start = as.integer(as.character(start)),
                        end = as.integer(as.character(end)),
                        gene = as.character(gene),
                        strand = as.character(strand) )
    
    return(de)
}
      

makeFastGenePos_GeneLevel <- function(file, outFile, keepXtrans = FALSE, keepNR = FALSE){

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

    saveRDS(rng, file = outFile)
    print(paste0("Saving gene info to ",outFile))
    return(rng)
}

