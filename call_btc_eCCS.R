#! /usr/bin/env Rscript

library(data.table)
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="", metavar="PATH"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="", metavar="PATH"),
  make_option(c("-m", "--model"), type="character", default=NULL, help="", metavar="PATH")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=TRUE)
opt <- parse_args(opt_parser)

file_path <- opt$input
model_path <- opt$model
output_path <- opt$output

source('./btc_CCS.fxn.R')

#file_path= "./data/tpm.txt"
#model_path="./data/LBR.tps.classifier.rds"
#output_path="./data/eCCS.txt"

if( all(file.exists(file_path) & file.exists(model_path)) ){
  
  classifier <- readRDS(model_path)
  rna_raw <- fread(file_path)
  
  rna_raw.mat <- as.matrix(rna_raw[,-c(1)])
  rownames(rna_raw.mat) <- rna_raw$gene_list
  
  predicted.w.score <- predict_TSP_with_confidence(tsp_classifier=classifier, tpm_matrix=rna_raw.mat)

  write.table( predicted.w.score[c('sample_id','rna_class','confidence.polarized')],  file = output_path, 
               row.names = F, col.names = T, 
               quote = FALSE, sep = "\t")

} else { cat("rna or model not found for ",file_path, "\n") 
  write.table( c(),  file = output_path, row.names = T, quote = FALSE, sep = "\t", col.names = F)
}
