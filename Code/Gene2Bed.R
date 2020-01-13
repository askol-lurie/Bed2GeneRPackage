#!/usr/bin Rscript

##############################
## Gene2Bed.R               ##
##                          ##
## Author: Andrew Skol      ##
##                          ##
## Created: 01/02/20        ##
## Updated: 01/13/20        ##
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
source(paste0(SCRDIR,'/Gene2Bed_funcs.R'))
##############

## GET GENE LOCATIONS USED TO CREATE RESOURCE LOCATION FILE ##
## geneLocsFile <- paste0(dirname(SCRDIR), "/Resources/GeneStartEnd37.rds")
geneLocsFile19.1 <- paste0(ResourceDir, "Genes_GenesPredictions_UCSCRefSeq_GRCh37.gz")
geneLocsFile19.2 <- paste0(ResourceDir, "Genes_GenesPredictions_NCBIRefSeq_GRCh37.gz")
geneLocsFile19.3 <- paste0(ResourceDir, "Genes_GenesPredictions_OtherUCSCRefSeq_GRCh37.gz")
geneLocsFile19.4 <- paste0(ResourceDir, "Genes_GenesPredictions_GENCODEV31lift37_Comprehensive_GRCh37.gz")
geneLocsFile19.5 <- paste0(ResourceDir, "Genes_GenesPredictions_GENCODEV31_psuedogenes_GRCh37.gz")
                           
geneLocsFile38.1 <- paste0(ResourceDir, "Genes_GenesPredictions_UCSCRefSeq_GRCh38.gz")
geneLocsFile38.2 <- paste0(ResourceDir, "Genes_GenesPredictions_NCBIRefSeq_GRCh38.gz")
geneLocsFile38.3 <- paste0(ResourceDir, "Genes_GenesPredictions_OtherUCSCRefSeq_GRCh38.gz")
geneLocsFile38.4 <- paste0(ResourceDir, "Genes_GenesPredictions_GENCODEV31_Comprehensive_GRCh38.gz")
geneLocsFile38.5 <- paste0(ResourceDir, "Genes_GenesPredictions_GENCODEV31_psuedogenes_GRCh38.gz")

geneLocsFiles19 <- c(geneLocsFile19.1, geneLocsFile19.2, geneLocsFile19.3, geneLocsFile19.4, geneLocsFile19.5)
geneLocsFiles38 <- c(geneLocsFile38.1, geneLocsFile38.2, geneLocsFile38.3, geneLocsFile38.4, geneLocsFile38.5)

## STORED GENE LOCATION FILE FOR FAST LOADING ##
geneLocFile19 <- paste0(ResourceDir,"GeneLocs_hg19.rds")
geneLocFile38 <- paste0(ResourceDir,"GeneLocs_hg38.rds")

## FILE WITH GENES NOT IN NCBI FILES ##
trickyGenesFile <- paste0(dirname(SCRDIR), "/Resources/TrickyGenesPositions.txt")

############

option_list = list(
  make_option(c("-b", "--build"), type="character", default=NULL, 
              help="Reference Build: hg19 or hg38",
              metavar="character"),
  make_option(c("-g", "--geneFile"), type="character", default=NULL,
              help="List of genes to assign bed intervals (gene in first column)",
              metavar="character"),
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
if (is.null(opt$build) ){
  print_help(opt_parser)
  stop("A reference build must be specified (hg19 or hg38)")
}else{
    if (opt$build %in% c("hg19","hg38") == FALSE){
        print_help(opt_parser)
        stop("Build must be either hg19 or hg38")
    }
}
## ENSURE GENE FILE(S) EXIST ##
if (!is.null(opt$geneFile)){
    FileMissing <- file.exists(opt$geneFile) == FALSE
    if( any(FileMissing)){
        stop(paste0("Input file : ",opt$geneFile[FileMissing],
                    " not found!"), call.=FALSE )
    }
}else{
    print_help(opt_parser)
    stop("Must specify a gene file.", call.=FALSE)
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



###################
## ONLY NEEDS TO BE RUN IF RECREATING GENE LOCATION FILES
###################
if (0){
    geneLocsFiles <- geneLocsFiles38
    if (opt$build == "hg19"){ geneLocsFiles <- geneLocsFiles19}
    makeGeneLocFile(files=geneLocsFiles, ResourceDir, Prefix = "GeneLocs", build = opt$build,
                    keepXtrans = FALSE, keepNR = TRUE)
}

if (0){
     geneLocsFiles <- geneLocsFiles38
    if (opt$build == "hg19"){ geneLocsFiles <- geneLocsFiles19}
    makeCodingLocFile(files=geneLocsFiles, ResourceDir, Prefix = "CodingLocs", build = opt$build,
                    keepXtrans = FALSE, keepNR = TRUE)
}
## GET BED DATA ##
genes <- getGenes(opt$geneFile)
outDir <- opt$out

## GET STORED GENE LOCATION INFORMATION
geneLocFile <- geneLocFile38
if (opt$build == "hg19"){
    geneLocFile <- geneLocFile19
}
geneLocs <- readRDS(geneLocFile)
geneLocs$gene <- as.character(geneLocs$gene)

prefix <- paste0(prefix,"_",opt$build)
gene2bed(genes, geneLocs, prefix, outDir)


             
    
    

