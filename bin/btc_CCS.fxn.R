#! /usr/bin/env Rscript

library(data.table)
library(optparse)

predict_TSP_with_confidence <- function(tsp_classifier, tpm_matrix, confidence_cutoff = 10/16) {
  
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
  predictions$rna_class[predictions$predicted_label == 0] <- 'eCCS-A'
  predictions$rna_class[predictions$predicted_label == 1] <- 'eCCS-B'
  
  predictions$confidence.polarized <- predictions$confidence
  predictions$confidence.polarized[predictions$predicted_label == 1 & !is.na(predictions$predicted_label)] <-
    1 - predictions$confidence[predictions$predicted_label == 1 & !is.na(predictions$predicted_label)]
  
  predictions$predicted_label[predictions$confidence <= confidence_cutoff] <- NA
  predictions$rna_class[predictions$confidence <= confidence_cutoff] <- NA
  predictions$sample_id <- colnames(tpm_matrix)
  
  return(predictions)
}


predict_gCCS <- function(cn, vaf, genebed, cytoband, full_model, THRESHOLD = 0.5 ){
  
  require(dplyr)
  require(tibble)
  require(CNTools)
  
  segmented.data <- CNSeg(cn)
  
  segment.gene <- getRS(
    segmented.data, 
    by="gene", 
    imput=FALSE, 
    XY=FALSE, 
    geneMap=genebed, 
    what="min")
  
  segment.gene <- rs(segment.gene)
  
  segment.gene.cyto <- inner_join(cytoband, segment.gene,  by=c('Hugo_Symbol'='genename'), relationship = "many-to-many")
  
  col.names <- names(segment.gene.cyto)[-c(1:7)]
  
  arm.med <- segment.gene.cyto %>%
    group_by(arm ) %>% 
    summarise(across(all_of(col.names), \(x) median(x, na.rm = TRUE)))
  
  transposed.raw <- arm.med %>%
    column_to_rownames(var = names(arm.med)[1]) %>%  
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sample") %>%
    mutate(across(-sample, as.numeric)) %>%  
    as_tibble()
  
  names(transposed.raw)[-1] <- paste0('arm.',names(transposed.raw)[-1])
  
  # add VAF
  transposed.raw <- inner_join(vaf, transposed.raw, by=c('sample'='sample'))
  
  # predict
  transposed.raw$prob <- predict(full_model, transposed.raw , type = "response")
  transposed.raw$glm_labels <- ifelse(transposed.raw$prob < THRESHOLD, 'gCCS-B', 'gCCS-A')
  transposed.raw$prob <- round(transposed.raw$prob, 3)
  transposed <- transposed.raw %>% dplyr::select(sample,glm_labels,prob)
  
  return(transposed)
}
