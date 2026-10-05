#! /usr/bin/env Rscript

library(data.table)
library(optparse)

check_columns <- function(df, required, what) {
  missing.cols <- setdiff(required, names(df))
  if (length(missing.cols) > 0) {
    stop('The ', what, ' is missing column(s): ', paste(missing.cols, collapse = ', '),
         '. Expected: ', paste(required, collapse = ', '))
  }
}

# Fail early with a usage message if any required option is absent or points at a missing file
check_input_files <- function(opt, required, parser) {
  absent <- required[vapply(required, function(x) is.null(opt[[x]]), logical(1))]
  if (length(absent) > 0) {
    print_help(parser)
    stop('Missing required option(s): ', paste0('--', absent, collapse = ', '), call. = FALSE)
  }
  not.found <- required[!file.exists(unlist(opt[required]))]
  if (length(not.found) > 0) {
    stop('File(s) not found: ',
         paste0('--', not.found, ' ', unlist(opt[not.found]), collapse = '; '), call. = FALSE)
  }
}

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
    
    for(this.gene in missing_genes){
      message(this.gene, ' is missing, trying aliases')

      # alias_symbol and prev_symbol are each '|'-separated lists (either may be empty or NA)
      these.aliases <- unlist(strsplit(c(all.aliases$alias_symbol[all.aliases$symbol == this.gene],
                                         all.aliases$prev_symbol[all.aliases$symbol == this.gene]), '\\|'))
      these.aliases <- these.aliases[!is.na(these.aliases) & these.aliases != '']

      found_gene <- intersect(these.aliases, rownames(tpm_matrix))
      if(length(found_gene) > 0){

        message(this.gene, ' found as ', found_gene[1])
        rownames(tpm_matrix)[rownames(tpm_matrix) == found_gene[1]] <- this.gene

      } else {
        message(paste0(this.gene," is missing from the TPM matrix and no aliases were found."))

        tsp_pairs <-
          tsp_pairs[!apply(tsp_pairs, 1, function(row) this.gene %in% row), , drop = FALSE]
        n_pairs <- nrow(tsp_pairs)
        warning(paste0(this.gene," has been removed. The algorithm was not trained this way. Use at own risk."))

      }
    }

    if (n_pairs == 0) {
      stop('None of the classifier gene pairs could be found in the expression table')
    }

  }

  predictions <- apply(tpm_matrix, 2, function(sample_expr) {

    votes <- sapply(1:n_pairs, function(i) {

      gene1 <- tsp_pairs[i, 1]
      gene2 <- tsp_pairs[i, 2]
      if (is.na(sample_expr[gene1]) || is.na(sample_expr[gene2])) {
        return(NA_character_)  # pair can't vote
      } else if (sample_expr[gene1] > sample_expr[gene2]) {
        return(labels[1])  # Class 0
      } else {
        return(labels[2])  # Class 1
      }

    })

    n_votes <- sum(!is.na(votes))
    if (n_votes == 0) {
      return(c(predicted_label = NA, confidence = NA))
    }

    vote_table <- table(votes)
    predicted_label <- names(which.max(vote_table))
    confidence <- max(vote_table) / n_votes
    
    # Set label to NA if confidence is exactly 0.5
    if (confidence == 0.5) {
      predicted_label <- NA
    }
    
    return(c(predicted_label = predicted_label, confidence = confidence))
  })
  
  n_na_votes <- colSums(is.na(tpm_matrix[unique(c(tsp_pairs)), , drop = FALSE]))
  if (any(n_na_votes > 0)) {
    warning('Missing (NA) expression for classifier genes in: ',
            paste(names(n_na_votes)[n_na_votes > 0], collapse = ', '),
            '. Those gene pairs were skipped for those samples.')
  }

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

  check_columns(cn, c('ID','chrom','loc.start','loc.end','num.mark','seg.mean'), 'segment file')
  check_columns(vaf, c('sample','mean.vaf'), 'VAF file')
  check_columns(genebed, c('chrom','start','end','geneid','genename'), 'gene bed')
  check_columns(cytoband, c('Hugo_Symbol','arm'), 'cytoband file')

  cn <- as.data.frame(cn)
  vaf <- as.data.frame(vaf)
  cn$ID <- as.character(cn$ID)
  vaf$sample <- as.character(vaf$sample)

  # segments and gene bed must use the same chromosome naming
  if (!any(grepl('^chr', cn$chrom)) && all(grepl('^chr', genebed$chrom))) {
    message('Adding "chr" prefix to segment chromosome names to match the gene bed')
    cn$chrom <- paste0('chr', cn$chrom)
  }

  if (sum(genebed$genename %in% cytoband$Hugo_Symbol) < 0.5 * length(unique(genebed$genename))) {
    stop('Fewer than half of the gene bed\'s genename values are in the cytoband file. ',
         'genename should hold HGNC gene symbols (check that the geneid/genename columns are not swapped).')
  }

  seg.samples <- unique(cn$ID)
  no.vaf <- setdiff(seg.samples, vaf$sample)
  if (length(no.vaf) == length(seg.samples)) {
    stop('No sample IDs are shared between the segment file (ID) and the VAF file (sample)')
  } else if (length(no.vaf) > 0) {
    warning('No mean VAF for: ', paste(no.vaf, collapse = ', '), '. These samples will not be classified.')
  }

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
  
  col.names <- intersect(as.character(unique(cn$ID)), names(segment.gene.cyto))
  
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
  model.terms <- attr(terms(full_model), 'term.labels')
  missing.terms <- setdiff(model.terms, names(transposed.raw))
  if (length(missing.terms) > 0) {
    stop('No genes from the segment file fall on: ', paste(missing.terms, collapse = ', '),
         '. Check that the segments and gene bed are on the same genome build.')
  }

  transposed.raw$prob <- predict(full_model, transposed.raw , type = "response")
  if (any(is.na(transposed.raw$prob))) {
    warning('Could not classify (missing arm-level copy number or VAF): ',
            paste(transposed.raw$sample[is.na(transposed.raw$prob)], collapse = ', '))
  }
  transposed.raw$glm_labels <- ifelse(transposed.raw$prob < THRESHOLD, 'gCCS-B', 'gCCS-A')
  transposed.raw$prob <- round(transposed.raw$prob, 3)
  transposed <- transposed.raw %>% dplyr::select(sample,glm_labels,prob)
  
  return(transposed)
}
