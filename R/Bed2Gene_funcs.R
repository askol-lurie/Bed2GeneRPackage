#' Annotate BED file with gene names
#'
#' Annotate BED file with gene names based on overlap between the BED interval and coding exons
#'
#' @param file BED file
#' @param genes Only output and annotation for intervals that overlap genes in list. Default is an empty vector which results in all genes being output. (optional)
#' @param geneFile File containing genes. Only output and annotation for intervals that overlap genes in the file will be output. Alternately, genes can be entered using the gene argument. (optional)
#' @param prefix Prefix for output files (optional)
#' @param outDir Output directory. If not given results will be printed to pwd(). (optional)
#' @param build Build reference for the position coordingates int he BED file (hg19 or hg38)
#' @param rmChrM Indicates if mitochondrial genes should be output and annotated (TRUE/FALSE)
#' 
#' @return Null
#'
#' @author Andrew Skol, \email{askol@@luriechildrens.org}
#' @seealso \code{\link{gene2bed}}
#' @keywords utilities
#'
#' @examples bed2gene(file = "MyBedFile.bed", prefilx = "myGenes.bed", rmChrM=TRUE)
#'
#' @export
#'        
bed2gene <- function(file, 
                     genes = c(), 
                     geneFile = "", 
                     prefix = "", 
                     outDir = pwd(), 
                     build = "hg19",
                     rmChrM=TRUE){

    checkArgs(file, genes, geneFile, prefix, outDir, build, rmChrM)
  
    ## SET UP OUTPUT FILE PREFIX ##
    if (prefix == ""){
        outFile <- basename(file)        
    }else{
        outFile <- prefix
    }
    outFile <- paste0(outFile , "_")

    ## set codinglocs based on build
    if (build == "hg19"){
      codingLocs <- exons19
      geneLocs <- genes19
      
    }else{
      codingLocs = exons38
      geneLocs <- genes38
    }
    
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

    ## determine if there is a header in the bed file (there shouldn't be)
    ## and read in accordingly
    tmp <- read.table(file = file, as.is=T, header=F, nrow=1, sep="\t")
    skipRow = 0
    if (any(tmp[,3] %in% c("chr","start","end","stop")) |
        any(sapply(tmp[1,2:3], class) %in% c("integer","numeric")) == FALSE){
        skipRow <- 1
    }
    
    bed <- read.table(file = file, as.is=T, header=F, sep="\t", skip=skipRow)
    names(bed)[1:3] <- c("chr","start","stop")
    bed$strand = "+"

    ## Remove any extra columns that were in the bed file ##
    bed <- bed %>% dplyr::select(chr, start, stop, strand) %>% 
      dplyr::mutate(chr = as.character(chr))
        
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
    bed <- bed %>% as.data.frame() 
    ## bed <- bed %>% dplyr::rename(chr = seqnames)
    
    FileOut <- paste0(outDir,"/", outFile, "AllInts_wGenes.bed")
    write.table(file = FileOut, bed, quote=F, row.names=F, col.names=T,
                sep="\t")
    print(paste0("Wrote bed file with gene names to ",FileOut))

    print(paste0("***** number of genes in genes ",length(genes)))
    if (length(genes) > 0){
        ## DETERMINE IF ANY GENES IN GENE LIST DO NOT HAVE BED INTERVALS ##
        missedGenesCoding <- genes[genes %in% unique(bed$gene_coding) == F]
        missedGenesTxn <- genes[genes %in% unique(bed$gene_txn) == F]
        
        if (length(missedGenesCoding) > 0){

            ## REPORT GENES IN GENE LIST WITH NO CODING REGIONS
            ## IN BED FILE
            missedGenesTbl <- as.data.frame(codingLocs[which(codingLocs$gene %in% missedGenesCoding)])
            
            FileOut <- paste0(outDir,"/", outFile, "GenesCodingWoIntervals.txt")
            write.table(file = FileOut, missedGenesTbl, quote=F, row.names=F, col.names=T, sep="\t")
            
            print(paste0(length(missedGenesCoding), " gene coding regions did not have intervals overlapping them!"))
            print(paste0(unique(missedGenesTbl$gene), collapse=", "))
            print(paste0("Missed genes coding regions are written to ",FileOut))

            ## REPORT GENES IN GENE LIST WITH REGIONS
            ## IN BED FILE (ANYWHERE IN GENE)
            missedGenesTbl <- as.data.frame(codingLocs[which(codingLocs$gene %in% missedGenesTxn)])
            
            FileOut <- paste0(outDir,"/", outFile, "GenesTxnWoIntervals.txt")
            write.table(file = FileOut, missedGenesTbl, quote=F, row.names=F, col.names=T, sep="\t")
            
            print(paste0(length(missedGenesTxn), " gene did not have intervals overlapping them!"))
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

checkArgs <-  function(file, genes, geneFile, prefix, outDir, build, rmChrM){
  
  if (is.null(file)){
    stop("You must specify a BED file (file = <bedfile.bed>", 
         call. = FALSE)
  }
  if (file.exists(file) == FALSE){
    stop(paste0("Your BED file was not found: ",file), call. = FALSE)
  }
  if (length(genes) != 0 & length(geneFile) != 0){
    stop("You must genes or geneFile, not both")
  }
  if (dir.exists(outDir) == FALSE){
    stop(paste("The specified output directory was not found:", outDir),
         call. = FALSE)
  }
  if (build %in% c("hg19", "hg38") == FALSE){
    stop("Build agruement must be hg19 or hg38")
  }
  if (as.character(rmChrM) %in% c("TRUE","FALSE","T","F") == FALSE){
    stop("rmChrM must be TRUE or FALSE", call. = FALSE)
  }
}


getGenes <- function(geneFiles){
  
  genes <- c()
  for (file in geneFiles){
    
    tmp <- read.table(file = file, as.is=T, header=FALSE)
    genes <- c(genes, c(unlist(tmp)))
    
  }
  return(genes)
}



