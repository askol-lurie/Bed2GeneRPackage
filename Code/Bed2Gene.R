#!/usr/bin Rscript

##############################
## Bed2Gene.R               ##
##                          ##
## Author: Andrew Skol      ##
##                          ##
## Created: 11/27/19        ##
## Updated: 11/27/19        ##
##                          ##
##############################
##
## DETERMINE DIRECTORY CODE IS LOCATED IN ##
## SUPPORTING *_FUNCS.r WILL BE IN SAME DIRECTORY ##
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
SCRDIR <- dirname(script.name)

## LOAD FUNCTIONS ##
source(paste0(SCRDIR,'/Bed2Gene_funcs.R'))
##############

## GET GENE LOCATIONS ##
geneLocsFile <- paste0(dirname(SCRDIR), "/Resources/GeneStartEnd37.rds")
geneLocs <- readRDS(geneLocsFile)
############

option_list = list(
  make_option(c("-b", "--bedFiles"), type="character", default=NULL, 
              help="Bed file or files (space seperated)" ,
              metavar="character"),
  make_option(c("-B", "--bedDir"), type="character", default=NULL, 
              help="Directory containing bed files to process",
              metavar="character"),
  make_option(c("-g", "--geneFiles"), type="character", default=NULL,
              help="List of genes to assign bed intervals to",
              metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="Output directory", metavar="character"))
  

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
    stop("Must specify one or more gene files.", call.=FALSE)
}
## ENSURE RUN DIRECTORY EXISTS 
if (!is.null(opt$PatDir)){
    if (!dir.exists(opt$bedDir)){
        stop(paste0("Input file : ",opt$bedDir, " not found!"), call.=FALSE )
    }
}
## ENSURE THE OUTPUT DIRCTORY EXISTS
if (dir.exists(dirname(opt$out)) == FALSE){
  stop(paste0("Output directory : ", dirname(opt$out), " not found!"), call.=FALSE )
}

## ################### ##
## END OPTION HANDLING ##
## ################### ##

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
    
    newBedFile <- bed2gene(file, genes, geneLocs, outDir)

}
             
    
    


    
