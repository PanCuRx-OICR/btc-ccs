#! /usr/bin/env Rscript

library(data.table)
library(optparse)

# locate the repository from this script's own path, so it can be run from any directory
script_arg <- grep('^--file=', commandArgs(trailingOnly = FALSE), value = TRUE)
repo_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub('^--file=', '', script_arg[1]))) else getwd()

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, metavar="PATH",
              help="Expression table (tab-separated): first column gene symbols, one column per sample. Trained on TPM. [required]"),
  make_option(c("-o", "--output"), type="character", default=NULL, metavar="PATH",
              help="Where to write the eCCS calls (tab-separated). [required]"),
  make_option(c("-m", "--model"), type="character", metavar="PATH",
              default=file.path(repo_dir, "data", "LBR.tps.classifier.rds"),
              help="Trained top-scoring-pair classifier. [default: %default]")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=TRUE,
                           description="Classify bulk tumour RNA into biliary tract cancer expression subtypes (eCCS-A / eCCS-B).")
opt <- parse_args(opt_parser)

source(file.path(repo_dir, 'bin', 'btc_CCS.fxn.R'))

check_input_files(opt, c('input', 'model'), opt_parser)
if (is.null(opt$output)) {
  print_help(opt_parser)
  stop('Missing required option: --output', call. = FALSE)
}

#opt$input= "./data/tpm.txt"
#opt$model="./data/LBR.tps.classifier.rds"
#opt$output="./data/eCCS.txt"

classifier <- readRDS(opt$model)
rna_raw <- fread(opt$input)

if (ncol(rna_raw) < 2) {
  stop('The expression table needs a gene column plus at least one sample column (is it tab-separated?)', call. = FALSE)
}

# first column holds gene symbols, whatever it is named
genes <- as.character(rna_raw[[1]])
rna_raw <- rna_raw[, -1, with = FALSE]

non.numeric <- names(rna_raw)[!vapply(rna_raw, is.numeric, logical(1))]
if (length(non.numeric) > 0) {
  stop('Non-numeric values in sample column(s): ', paste(non.numeric, collapse = ', '), call. = FALSE)
}

dup.genes <- intersect(unique(genes[duplicated(genes)]), c(classifier$TSPs))
if (length(dup.genes) > 0) {
  warning('Classifier gene(s) appear more than once; using the first row for: ', paste(dup.genes, collapse = ', '))
}

rna_raw.mat <- as.matrix(rna_raw)
rownames(rna_raw.mat) <- genes
rna_raw.mat <- rna_raw.mat[!duplicated(genes), , drop = FALSE]

predicted.w.score <- predict_TSP_with_confidence(tsp_classifier=classifier, tpm_matrix=rna_raw.mat)

write.table( predicted.w.score[c('sample_id','rna_class','confidence.polarized')],  file = opt$output,
             row.names = F, col.names = T,
             quote = FALSE, sep = "\t")
