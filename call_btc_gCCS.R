#! /usr/bin/env Rscript

library(data.table)
library(optparse)

option_list = list(
  make_option(c("-s", "--seg_path"), type="character", default=NULL, help="", metavar="PATH"),
  make_option(c("-v", "--vaf_path"), type="character", default=NULL, help="", metavar="PATH"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="", metavar="PATH"),
  make_option(c("-m", "--model"), type="character", default=NULL, help="", metavar="PATH"),
  make_option(c("-b", "--bed_path"), type="character", default=NULL, help="", metavar="PATH"),
  make_option(c("-c", "--cytoband_path"), type="character", default=NULL, help="", metavar="PATH")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=TRUE)
opt <- parse_args(opt_parser)

seg_path <- opt$seg_path
model_path <- opt$model
output_path <- opt$output
vaf_path <- opt$vaf_path
bed_path <- opt$bed_path
cytoband_path <- opt$cytoband_path

source('./bin/btc_CCS.fxn.R')

#seg_path= "./data/example.seg"
#vaf_path= "./data/mean_vaf.txt"
#model_path="./data/cnv.model.rds"
#cytoband_path="./data/cytoband.txt"
#genome_build='hg38'
#bed_path=paste0("./data/genebed.",genome_build,".txt")
#output_path="./data/gCCS.txt"

if( all(file.exists(seg_path) & file.exists(vaf_path) &
        file.exists(model_path) & file.exists(cytoband_path) &
        file.exists(bed_path)) ){
  
  full_model <- readRDS(model_path)
  vaf <- fread(vaf_path)
  cn <- fread(seg_path)
  cytoband <- fread(cytoband_path)
  genebed <- fread(bed_path)
  
  gCCS_predictions <- predict_gCCS(cn = cn, 
                                   vaf = vaf, 
                                     genebed, cytoband, full_model )
  
  
  write.table( gCCS_predictions,  file = output_path, 
               row.names = F, col.names = T, 
               quote = FALSE, sep = "\t" )

} else { 
  message("some files not found") 
  write.table( c(),  file = output_path, row.names = T, quote = FALSE, sep = "\t", col.names = F)
}
