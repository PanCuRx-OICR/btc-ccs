#! /usr/bin/env Rscript

library(data.table)
library(optparse)

# locate the repository from this script's own path, so it can be run from any directory
script_arg <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
repo_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub('^--file=', '', script_arg[1]))) else getwd()

option_list = list(
  make_option(c("-s", "--seg_path"), type="character", default=NULL, metavar="PATH",
              help="Copy-number segments (tab-separated): ID, chrom, loc.start, loc.end, num.mark, seg.mean. [required]"),
  make_option(c("-v", "--vaf_path"), type="character", default=NULL, metavar="PATH",
              help="Mean variant allele frequency per sample (tab-separated): sample, mean.vaf. [required]"),
  make_option(c("-o", "--output"), type="character", default=NULL, metavar="PATH",
              help="Where to write the gCCS calls (tab-separated). [required]"),
  make_option(c("-m", "--model"), type="character", metavar="PATH",
              default=file.path(repo_dir, "data", "cnv.model.rds"),
              help="Trained copy-number model. [default: %default]"),
  make_option(c("-b", "--bed_path"), type="character", metavar="PATH",
              default=file.path(repo_dir, "data", "genebed.hg38.txt"),
              help="Gene coordinates on the same build as the segments (genebed.hg38.txt or genebed.hg19.txt). [default: %default]"),
  make_option(c("-c", "--cytoband_path"), type="character", metavar="PATH",
              default=file.path(repo_dir, "data", "cytoband.txt"),
              help="Gene to chromosome-arm map. [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=TRUE,
                           description="Classify tumour copy number into biliary tract cancer genomic subtypes (gCCS-A / gCCS-B).")
opt <- parse_args(opt_parser)

source(file.path(repo_dir, 'bin', 'btc_CCS.fxn.R'))

check_input_files(opt, c('seg_path', 'vaf_path', 'model', 'bed_path', 'cytoband_path'), opt_parser)
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop('Missing required option: --output', call. = FALSE)
}

#opt$seg_path= "./data/example.seg"
#opt$vaf_path= "./data/mean_vaf.txt"
#opt$model="./data/cnv.model.rds"
#opt$cytoband_path="./data/cytoband.txt"
#opt$bed_path="./data/genebed.hg38.txt"
#opt$output="./data/gCCS.txt"

full_model <- readRDS(opt$model)
vaf <- fread(opt$vaf_path)
cn <- fread(opt$seg_path)
cytoband <- fread(opt$cytoband_path)
genebed <- fread(opt$bed_path)

gCCS_predictions <- predict_gCCS(cn = cn,
                                 vaf = vaf,
                                 genebed, cytoband, full_model )

write.table( gCCS_predictions,  file = opt$output,
             row.names = F, col.names = T,
             quote = FALSE, sep = "\t" )
