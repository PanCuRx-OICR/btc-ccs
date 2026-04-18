#! /usr/bin/env Rscript

library(data.table)
library(optparse)

predict_TSP_with_confidence <- function(tsp_classifier, tpm_matrix, confidence_cutoff = 0.8) {
  
  if(length(tsp_classifier$TSPs) > nrow(tsp_classifier$aliases['symbol'])){
    stop('There is a formatting error in aliases, some genes are missing')
  }
  
  tsp_pairs <- tsp_classifier$TSPs
  labels <- tsp_classifier$labels
  n_pairs <- nrow(tsp_pairs)
  
  # Ensure gene names are in the TPM matrix
  all_genes <- unique(c(tsp_pairs[,1], tsp_pairs[,2]))
  missing_genes <- setdiff(all_genes, rownames(tpm_matrix))
  
  if (length(missing_genes) > 0) {
    require(stringr)
    
    aliases <- tsp_classifier[['aliases']]
    all.aliases <- aliases[aliases$symbol %in% missing_genes, ]
    
    for(this.gene in all.aliases$symbol){
      message(this.gene, ' is missing, trying aliases')
      
      these.aliases <- paste0(all.aliases$alias_symbol[all.aliases$symbol == this.gene],
                              all.aliases$prev_symbol[all.aliases$symbol == this.gene], sep="|")
      
      these.aliases <- strsplit(these.aliases, '\\|')[[1]]
      found_gene <- intersect(these.aliases, rownames(tpm_matrix))
      if(length(found_gene) > 0){
        
        rownames(tpm_matrix)[rownames(tpm_matrix) == found_gene] <- this.gene
        
      } else {
        message(paste0(this.gene," is missing from the TPM matrix and no aliases were found."))
        
        tsp_pairs <- 
          tsp_pairs[!apply(tsp_pairs, 1, function(row) this.gene %in% row), ]
        n_pairs <- nrow(tsp_pairs)
        warning(paste0(this.gene," has been removed. The algorithm was not trained this way. Use at own risk."))
        
      }   
    }
    
    
  }
  
  predictions <- apply(tpm_matrix, 2, function(sample_expr) {
    votes <- sapply(1:n_pairs, function(i) {
      gene1 <- tsp_pairs[i, 1]
      gene2 <- tsp_pairs[i, 2]
      if (sample_expr[gene1] > sample_expr[gene2]) {
        return(labels[1])  # Class 0
      } else {
        return(labels[2])  # Class 1
      }
    })
    
    vote_table <- table(votes)
    predicted_label <- names(which.max(vote_table))
    confidence <- max(vote_table) / n_pairs
    
    
    # Set label to NA if confidence is exactly 0.5
    if (confidence == 0.5) {
      predicted_label <- NA
    }
    
    return(c(predicted_label = predicted_label, confidence = confidence))
  })
  
  # Transpose to get samples as rows
  predictions <- t(predictions)
  predictions <- as.data.frame(predictions)
  predictions$confidence <- as.numeric(predictions$confidence)
  
  predictions$rna_class <- NA
  predictions$rna_class[predictions$predicted_label == 0] <- 'CMS-A'
  predictions$rna_class[predictions$predicted_label == 1] <- 'CMS-B'
  
  predictions$confidence.polarized <- predictions$confidence
  predictions$confidence.polarized[predictions$predicted_label == 1 & !is.na(predictions$predicted_label)] <-
    1 - predictions$confidence[predictions$predicted_label == 1 & !is.na(predictions$predicted_label)]
  
  predictions$predicted_label[predictions$confidence <= confidence_cutoff] <- NA
  predictions$rna_class[predictions$confidence <= confidence_cutoff] <- NA
  
  return(predictions)
}

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

#file_path= paste0("/Volumes/pcsi/users/fbeaudry/tpm.txt")
#model_path=paste0("/Volumes/pcsi/references/cohort_lists/LBR.tps.classifier.rds")
#
if( all(file.exists(file_path) & file.exists(cohort_path) & file.exists(gene_list_path)) ){
  
  classifier <- readRDS(model_path)
  rna_raw <- fread(file_path)
  
  rna_raw.mat <- as.matrix(rna_raw[,-c(1,3)])
  rownames(rna_raw.mat) <- rna_raw$gene_name
  
  predicted.w.score <- predict_TSP_with_confidence(tsp_classifier=classifier, tpm_matrix=rna_raw.mat)
  write.table( predicted.w.score,  file = output_path, row.names = F, quote = FALSE, sep = "\t", col.names = T)

} else { cat("rna not found for ",file_path, "\n") 
  write.table( c(),  file = output_path, row.names = T, quote = FALSE, sep = "\t", col.names = F)
}
