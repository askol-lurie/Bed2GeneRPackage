
exons19 <- readRDS("GeneExons_hg19.rds")
exons38 <- readRDS("GeneExons_hg38.rds")
genes19 <- readRDS("GeneLocs_hg19.rds")
genes38 <- readRDS("GeneLocs_hg38.rds")

usethis::use_data(exons19, exons38, genes19, genes38, internal = TRUE)
