
library(optparse)
library(GenomicRanges)
library(dplyr)
library(data.table)

getGenes <- function(geneFiles){

    genes <- c()
    for (file in geneFiles){
        
        tmp <- read.table(file = file, as.is=T, header=FALSE)
        genes <- c(genes, c(unlist(tmp)))
        
    }
    return(genes)
}
        
bed2gene <- function(file, genes = c(), geneLocs, prefix = "", outDir){

    if (length(genes) > 0){
        
        ## FIRST SEE IF THERE ARE ANY MISSING GENES IN GENELOCS ##
        miss <- genes[genes %in% geneLocs$gene == FALSE]
        if (length(miss) > 0){
            
            fileMiss <- paste0(outDir, "/MissGenes.txt")
            print(paste0("The following genes were not found in the reference gene position list :"))
            print(miss)
            print(paste0("List of unfound genes printed to ", fileMiss))
            write.table(file = fileMiss, miss, quote=F, row.names = F, col.names=F)
        }
        
        ## SEARCH ONLY FOR THE GENES IN GENES ##
        geneLocs <- geneLocs[which(geneLocs$gene %in% genes)]
    }

    ## GET BED AND CONVERT TO GRANGE OBJECT
    bed <- read.table(file = file, as.is=T, header=T)
    names(bed)[1:3] <- c("chr","start","stop")
    bed$strand = "+"
    
    bed <- makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE, seqnames.field = "chr", 
                                    start.field = "start", end.field = "stop",
                                    strand.field = "strand")
    
    ## FIND THE GENE(S) THAT AN INTERVAL OVERLAPS WITH
    ol <- as.data.frame(findOverlaps(bed,geneLocs, type="any", ignore.strand=TRUE ))

    ol <- ol %>% mutate(gene = geneLocs$gene[subjectHits])

    ## remove multiple lines for each unique inverval-gene entry
    ol <- ol %>% group_by(queryHits, gene) %>% filter(row_number() == 1)

    ## if multiple genes map to the same interval then combined gene names and
    ## remove multiple occurance of interval
    ol <- ol %>% group_by(queryHits) %>% mutate( ifelse(n() > 1, paste(gene,collapse=";"), gene) )%>%
        filter(row_number() == 1)
        
    ## ASSIGN GENE(S) TO INTERVALS
    bed$gene_pred = ""
    bed$gene_pred[ol$queryHits] <- ol$gene

    if (prefix != ""){
        outFile <- basename(file)        
        outFile <- gsub("\\..+","", outFile)
    }else{
        outFile <- prefix
    }
    outFile <- paste0(outFile , "_")
    
    file <- paste0(outDir,"/", outFile, "AllInts_wGenes.bed")
    write.table(file = file, bed, quote=F, row.names=F, col.names=F,
                sep="\t")
    print(paste0("Wrote bed file with gene names to ",file))
    
    if (length(genes) > 0){
        ## DETERMINE IF ANY GENES IN GENE LIST DO NOT HAVE BED INTERVALS ##
        missedGenes <- genes[genes %in% unique(bed$gene) == F]
        
        if (length(missedGenes) > 0){
            
            missedGenesTbl <- as.data.frame(geneLocs[which(geneLocs$gene %in% missedGenes)])
            
            file <- paste0(outDir,"/", outFile, "GenesWoIntervals.txt")
            write.table(file = file, missedGenesTbl, quote=F, row.names=F, col.names=T, sep="\t")
            
            print(paste0(length(missedGenes), " gene did not have intervals overlapping them!"))
            print(paste0(unique(missedGenesTbl$gene), collapse=", "))
            print(paste0("Missed genes are written to ",file))
        }                
    }

    ## IF BED FILE CONTAINED GENES (COLUMN NAME gene)
    ## then also output bed intervals with those genes (given name) only ##
    bed <- bed %>% as.data.frame()
    if (any(names(bed) == "gene")){

        bed <- bed %>% rename(gene_given = gene) %>% filter(gene_given %in% genes)
        file <- paste0(outDir,"/", outFile, "GivenGenesOnly.bed")
        write.table(file = file, bed, quote=F, row.names=F, col.names=T, sep="\t")
        print("Writing bed intervals containing genes from gene file")
        print(paste0("Wrote ",file))
    }    
   
    return()
}

makeFastGenePos <- function(file, outFile, keepXtrans = FALSE, keepNR = FALSE){

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
    de <- exonify(d)
    
    rng <- makeGRangesFromDataFrame(de, keep.extra.columns = T, seqnames.field = "chr",
                                    start.field = "start", end.field = "end",
                                    strand.field = "strand")
    saveRDS(rng, file = outFile)
    print(paste0("Saving gene info to ",outFile))
}

exonify <- function(data){

    de <- list()
    for (i in 1:nrow(data)){
        if (i%%5000 == 0){
            print(paste0("Working on transcript ",i, " of ",nrow(data)))
        }
        tran = data$name[i]
        gene = data$name2[i]
        strand = data$strand[i]
        chrom = data$chrom[i]
        tmp <- cbind(strsplit(data$exonStarts[i], split=",")[[1]],
                     strsplit(data$exonEnds[i], split=",")[[1]])
        de[[i]] <- cbind(chrom, tmp, gene, strand, tran)
    }
    
    de <- do.call(rbind, de)

    de <- as.data.frame(de, string.as.factors = FALSE)    
    names(de) <- c("chr","start","end","gene","strand","tran")

    ## keep distinct intervals (means will delete transcript names) ##
    print("Removing duplicate intervals (as a result of mult transcripts per gene)")
    de <- de %>% group_by(chr,start,end,gene) %>% filter(row_number() == 1) %>% ungroup()
    
    return(de)
}
      

if (0){

    ## This is how the gene position file was made ##
    
    ## RESOURCE FILES ##
    ## ResourceDir <- "/home/askol/Projects/Bed2Gene/Resources/"
    ## UCSCGene37fast <- paste0(ResourceDir,"GeneStartEnd37.rds")
    UCSCGene37fast <- paste0(ResourceDir,"GeneExons37.rds")
    UCSCGeneFile37 <- paste0(ResourceDir,"Genes_GenesPredictions_UCSCRefSeq_GRCh37.gz")

    Gene37 <- makeFastGenePos(UCSCGeneFile37, UCSCGene37fast, keepXtrans = FALSE)
}
