#' Create a BED file containing the coding regions of genes
#'
#' Create a BED file of coding exons for a user-specified set of genes. Padding can be added#'
#' @param genes List of genes for which coding exons should be output coding exons. 
#' @param genefile A file that contains one gene per line
#' @param prefix Prefix for output files (optional)
#' @param outDir Output directory (optional)
#' @param pad Amount of padding up and downstream of coding exon (in bp) (default = 0)
#' @param rmChrM Indicates if mitochondrial genes should be output (TRUE/FALSE) (default is TRUE)
#' 
#' @return Null
#'
#' @author Andrew Skol, \email{askol@@luriechildrens.org}
#' @seealso \code{\link{bed2gene}}
#' @keywords utilities
#'
#' @examples gene2bed(c("TP53","HNF4A"), prefix="MyGenes", pad=10, rmChrM=TRUE)
#'
#' @export
#'            
gene2bed <- function(genes, geneFile, prefix, outDir, pad=0, rmChrM=TRUE){
  
  checkArgs_gene2bed(file, genes, geneFile, prefix, outDir, pad, rmChrM)
  
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

checkArgs_gene2bed <-  function(file, genes, geneFile, prefix, outDir, pad, rmChrM){
  
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
  if (as.character(rmChrM) %in% c("TRUE","FALSE","T","F") == FALSE){
    stop("rmChrM must be TRUE or FALSE", call. = FALSE)
  }
}






