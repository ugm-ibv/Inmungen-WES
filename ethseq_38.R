##Requirements
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

#BiocManager::install("SNPRelate", lib="/home/jperez/R")
#BiocManager::install("gdsfmt", lib="/home/jperez/R")
#install.packages("EthSEQ", repos = "https://cloud.r-project.org/")

##EthSEQ
library(gdsfmt, lib.loc="/home/jperez/R")
library(EthSEQ, lib.loc="/home/jperez/R")
library(SNPRelate, lib.loc="/home/jperez/R")

dir<-Sys.getenv("outDir")

setwd(paste0(dir, "/ethseq"))
download.file("https://github.com/cibiobcg/EthSEQ_Data/blob/master/EthSEQ_Models_hg38/Exonic.All.Model.hg38.gds", destfile = "/home/jperez/COVID19/results/pipeline/ethseq/Exonic.All.Model.hg38.gds")

## Construir reference model
#ethseq.RM(vcf.fn, annotations, out.dir = "./",
#          model.name = "Reference.Model", bed.fn = NA, call.rate = 1, cores = 1)

## Executar amb Exonic.all com a referència
ethseq.Analysis(
            target.vcf = "recalibrated_variants_DP20_gatk_het_ALL_vep_modified.vcf",
            model.gds = "Exonic.All.Model.hg38.gds",
            verbose=TRUE,
            out.dir=paste0(dir, "/ethseq"),
            composite.model.call.rate = 0.99,
            space = "3D", 
            cores = 4)

