#!/usr/bin Rscript

##############################
## Bed2Gene.R               ##
##                          ##
## Author: Andrew Skol      ##
##                          ##
## Created: 11/27/19        ##
## Updated: 01/15/20        ##
##                          ##
##############################
##
## DETERMINE DIRECTORY CODE IS LOCATED IN ##
## SUPPORTING *_FUNCS.r WILL BE IN SAME DIRECTORY ##
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
SCRDIR <- dirname(script.name)
ResourceDir <- paste0(dirname(SCRDIR),"/Resources/")

## LOAD FUNCTIONS ##
source(paste0(SCRDIR,'/Bed2Gene_funcs.R'))
##############

## GET GENE LOCATIONS USED TO CREATE EXON LOCATION FILES ##
geneLocsFile19.2 <- paste0(ResourceDir, "Genes_GenesPredictions_UCSCRefSeq_GRCh37.gz")
geneLocsFile19.1 <- paste0(ResourceDir, "Genes_GenesPredictions_NCBIRefSeq_GRCh37.gz")
geneLocsFile19.3 <- paste0(ResourceDir, "Genes_GenesPredictions_OtherUCSCRefSeq_GRCh37.gz")
geneLocsFile19.4 <- paste0(ResourceDir, "Genes_GenesPredictions_GENCODEV31lift37_Comprehensive_GRCh37.gz")
geneLocsFile19.5 <- paste0(ResourceDir, "Genes_GenesPredictions_GENCODEV31_psuedogenes_GRCh37.gz")

geneLocsFile38.2 <- paste0(ResourceDir, "Genes_GenesPredictions_UCSCRefSeq_GRCh38.gz")
geneLocsFile38.1 <- paste0(ResourceDir, "Genes_GenesPredictions_NCBIRefSeq_GRCh38.gz")
geneLocsFile38.3 <- paste0(ResourceDir, "Genes_GenesPredictions_OtherUCSCRefSeq_GRCh38.gz")
geneLocsFile38.4 <- paste0(ResourceDir, "Genes_GenesPredictions_GENCODEV31_Comprehensive_GRCh38.gz")
geneLocsFile38.5 <- paste0(ResourceDir, "Genes_GenesPredictions_GENCODEV31_psuedogenes_GRCh38.gz")
geneLocsFile38.6 <- paste0(ResourceDir, "Mitochondrial_genes.bed") ## copied from /home/win.ngs/NGS/medex_illumina/Bed_files

mitoFile <- paste0(ResourceDir, "Mitochondrial_genes.bed") ## copied from /home/win.ngs/NGS/medex_illumina/Bed_files
geneLocsFiles19 <- c(geneLocsFile19.1, geneLocsFile19.2, geneLocsFile19.3, geneLocsFile19.4, geneLocsFile19.5)
geneLocsFiles38 <- c(geneLocsFile38.1, geneLocsFile38.2, geneLocsFile38.3, geneLocsFile38.4, geneLocsFile38.5)

## GET GENE LOCATIONS ##
## geneLocsFile <- paste0(dirname(SCRDIR), "/Resources/GeneStartEnd37.rds")
exonLocsFile.19 <- paste0(ResourceDir, "GeneExons_hg19.rds")
exonLocsFile.38  <- paste0(ResourceDir, "GeneExons_hg38.rds")



############

option_list = list(
  make_option(c("-b", "--bedFiles"), type="character", default=NULL, 
              help="Bed file or files (space seperated) (-b or -B required)" ,
              metavar="character"),
  make_option(c("-B", "--bedDir"), type="character", default=NULL, 
              help="Directory containing bed files to process (-b or -B required)",
              metavar="character"),
  make_option(c("-g", "--geneFiles"), type="character", default=NULL,
              help="List of genes to assign bed intervals to (optional)",
              metavar="character"),
  make_option(c("-G", "--genome"), type="character", default=NULL, 
              help="Reference Genome Build: hg19 or hg38", metavar = "character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Output directory", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="Output File Prefix (optional)", metavar="character"))
  

## ################# ##
## OPTION HANDLING   ##
## ################# ##

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## ENSURE ARGUMENTS ARE GIVEN
if (is.null(opt$out)){
  print_help(opt_parser)
  stop("Input and output files are required", call.=FALSE)
}
## ENSURE NOT BOTH FILE AND DIRECTORY ARE SPECIFIED
if (!is.null(opt$bedFiles) & !is.null(opt$bedDir)){
  print_help(opt_parser)
  stop("Must either bed file(s) or a bed directory, but not both")
}
## ENSURE THE INPUT FILE EXISTS
if (!is.null(opt$bedFiles)){
    FileMissing <- file.exists(opt$bedFiles) == FALSE
    if( any(FileMissing)){
        stop(paste0("Input file : ",opt$bedFiles[FileMissing],
                    " not found!"), call.=FALSE )
    }
}
## ENSURE GENE FILE(S) EXIST ##
if (!is.null(opt$geneFiles)){
    FileMissing <- file.exists(opt$geneFiles) == FALSE
    if( any(FileMissing)){
        stop(paste0("Input file : ",opt$geneFiles[FileMissing],
                    " not found!"), call.=FALSE )
    }   
}else{
    print("No gene file provided. Examining all possible genes intervals may fall on.")
    ## stop("Must specify one or more gene files.", call.=FALSE)
}
## ENSURE GENOME BUILD IS SPECIFIED
if (!is.null(opt$genome)){
    if (opt$genome %in% c("hg19","hg38") == F){
        stop("Genome build must be either hg19 or hg38")
    }
}else{
    stop("A genome build must be supplied (hg19 or hg38)")
}
## ENSURE RUN DIRECTORY EXISTS 
if (!is.null(opt$PatDir)){
    if (!dir.exists(opt$bedDir)){
        stop(paste0("Input file : ",opt$bedDir, " not found!"), call.=FALSE )
    }
}
prefix = ""
if (!is.null(opt$prefix)){
    prefix = opt$prefix
}

## ENSURE THE OUTPUT DIRCTORY EXISTS
if (dir.exists(dirname(opt$out)) == FALSE){
  stop(paste0("Output directory : ", dirname(opt$out), " not found!"), call.=FALSE )
}

## ################### ##
## END OPTION HANDLING ##
## ################### ##

if (0){
    geneLocsFiles <- geneLocsFiles38
    if (opt$genome == "hg19"){ geneLocsFiles <- geneLocsFiles19}
    makeExonLocFile(files=geneLocsFiles, mitoFile = mitoFile, ResourceDir, Prefix = "GeneExons",
                    build = opt$genome,  keepXtrans = FALSE, keepNR = TRUE)
}

exonLocsFile <- exonLocsFile.38
if (opt$genome == "hg19"){
    exonLocsFile <- exonLocsFile.19
}
exonLocs <- readRDS(exonLocsFile)
exonLocs$gene <- as.character(exonLocs$gene)

## GET BED DATA ##
bedFiles <- geneFiles <- c()
if (!is.null(opt$bedFile)){
    bedFiles <- opt$bedFile
}else{
    bedFiles <- dir(opt$bedDir)
    if (length(bedFiles) == 0){
        stop(paste0("No files found in ",opt$bedDir), call.=FALSE )
    }        
}

genes <- getGenes(opt$geneFiles)
outDir <- opt$out

for (file in bedFiles){
    
    bed2gene(file, genes, exonLocs, prefix, outDir)

}
             
    
    


    
