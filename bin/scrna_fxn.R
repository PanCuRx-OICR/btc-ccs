
convert_w_to_tpm_like <- function(w.matrix, gene_set){
  
  gene_set$kb <- gene_set$Length / 1000
  
  gene_kb <- gene_set$kb
  names(gene_kb) <- gene_set$gene_name
  
  #rownames(w.matrix) <- w.matrix$ensid
  
  
  W_matrix <- w.matrix %>% 
    data.frame() %>%
    dplyr::select(starts_with("sig")) %>%
    tibble::rownames_to_column(var = "gene_name") %>%
    inner_join(., gene_set, by = c('gene_name'='gene_name')) %>%
    mutate(Gene_Symbol = make.unique(gene_name)) %>%
    column_to_rownames(var = "Gene_Symbol") %>% 
    dplyr::select(-all_of(c('gene_name','Chr','Start','End','strand','biotype', 'Length', 'ensid','kb')))
  
  W_matrix <- W_matrix[complete.cases(W_matrix),]
  
  tpm_data <- exp(W_matrix) - 1
  
  tpm_data <- tpm_data / gene_kb[rownames(W_matrix)]
  scaling_factors <- colSums(tpm_data, na.rm = T)
  
  tpm_data <- (tpm_data/ scaling_factors) *  1e6
  
  tpm_data <- log(tpm_data + 1)
  tpm_data <- tpm_data[complete.cases(tpm_data),]
  
  return(tpm_data)
}

get_cms_from_pseudobulk <- function(tumour.subset_obj, confidence_cutoff  = 0.6, 
                                    grouping_var = "sample_id", 
                                    basedir = '~/Documents/scripts/github/btc-ccs/'){
  
  source(file.path(basedir,"/bin/BTC.functions.R"))
  classifier <- readRDS(file.path(basedir,'/results/tps.classifier.rds'))
  
  counts <- GetAssayData(tumour.subset_obj, slot = "counts")
  
  cell_metadata <- tumour.subset_obj@meta.data
  groupings <- split(colnames(tumour.subset_obj), cell_metadata[[grouping_var]])
  
  # Create a list of summed counts per group
  pseudobulk_list <- lapply(groupings, function(cells) {
    Matrix::rowSums(counts[, cells, drop = FALSE])
  })
  
  # Combine into a matrix
  pseudobulk_matrix <- do.call(cbind, pseudobulk_list)
  colnames(pseudobulk_matrix) <- names(pseudobulk_list)
  
  pseudobulk_df <- data.frame(pseudobulk_matrix)
  pseudobulk_df$gene_name <- row.names(pseudobulk_matrix)
  
  pseudobulk_df <- left_join(pseudobulk_df, gene_set, by=c('gene_name'='gene_name')) %>% 
    filter(!is.na(ensid))
  
  n.samples <- tumour.subset_obj@meta.data %>% dplyr::select(all_of(grouping_var)) %>% 
    unique() %>% nrow()
  
  pseudo.rna.mat <- apply(as.matrix(pseudobulk_df[,c(1:n.samples)] ), 2, as.numeric)
  rpk  <- (pseudo.rna.mat / (pseudobulk_df$Length / 1000)) 
  scaling_factors <- colSums(rpk)
  tpm <- t(t(rpk) / scaling_factors) * 1e6
  rownames(tpm) <- pseudobulk_df$gene_name
  
  predicted.w.score <- 
    predict_TSP_with_confidence(tpm_matrix=tpm, 
                                tsp_classifier=classifier, 
                                confidence_cutoff = confidence_cutoff )
  
  return(predicted.w.score)
  
}

get_cell_marker_list <- function(basedir = '~/Documents/scripts/github/btc-ccs/'){
  # Grab cell type markers from Cellmarker 2.0 
  
  cell_marker_db <- file.path(basedir,'/source_data/Cell_marker_Seq.xlsx')
  
  cell_markers <- read.xlsx(cell_marker_db) %>%
    data.frame() %>%
    dplyr::select(species, cancer_type,
                  cell_type, tissue_type,
                  cell_name,marker, Symbol,
                  Genetype, Genename) %>%
    dplyr::filter(species == "Human") %>%
    dplyr::filter(!is.na(Symbol), Symbol != "", 
                  !is.na(cell_name), cell_name != "") %>%
    dplyr::distinct(cell_name, Symbol) 
  
  # convert cell marker df into a list
  
  cell_marker_list <- cell_markers %>%
    split(.$cell_name) %>%
    lapply(function(df) unique(df$Symbol))
  
  return(cell_marker_list)
}

get_overlap_scores <- function(query_genes, reference_list) { 
  
  total_genes <- unique(unlist(reference_list)) %>% union(query_genes)
  
  result <- lapply(reference_list, function(ref_genes) {
    
    a <- length(intersect(query_genes, ref_genes))                    # overlap
    b <- length(setdiff(query_genes, ref_genes))                      # only in query
    c <- length(setdiff(ref_genes, query_genes))                      # only in ref
    d <- length(setdiff(total_genes, union(query_genes, ref_genes)))  # in neither
    mat <- matrix(c(a, b, c, d), nrow = 2)
    pval <- fisher.test(mat, alternative = "greater")$p.value
    return(pval)
  })
  
  result_df <- data.frame(cell_type = names(result), pval = unlist(result))
  result_df <- result_df %>% arrange(pval)
  return(result_df)
}

make_gene_sets <- function(tpm_data, ntile=0.8){
  suppressPackageStartupMessages({
    require(mclust)
    require(lsa)
  })
  
  gene_sets <- list()
  for(this.sig.index in c(1:ncol(tpm_data))){
    
    cat(paste0('sig.',this.sig.index),'\n')
    
    gene_mat <- as.matrix(tpm_data[order(unlist(tpm_data[paste0('sig.',this.sig.index)]), decreasing = T),][paste0('sig.',this.sig.index)])
    gene_list <- as.double(gene_mat)
    names(gene_list) <- rownames(gene_mat)
    
    gene.df <- data.frame(gene_list)
    clustered.weights <- Mclust(log(gene.df$gene_list), G=2)
    gene.df$class <- clustered.weights$classification
    gene.df$gene <- row.names(gene.df)
    #print(ggplot(gene.df, aes(x=gene_list, fill=as.factor(class))) + geom_histogram() + theme_bw())
    
    
    gene.df <- gene.df[gene.df$class == 2,]
    
    top_quartile_threshold <- quantile(gene.df$gene_list, ntile)
    
    genelist.filt <- gene.df$gene_list[gene.df$gene_list >= top_quartile_threshold]
    names(genelist.filt) <- gene.df$gene[gene.df$gene_list >= top_quartile_threshold]
    gene_sets[[paste0('sig.',this.sig.index)]] <-  names(genelist.filt)
    
  }
  
  
  # Create a binary presence matrix
  all_genes <- unique(unlist(gene_sets))
  binary_matrix <- matrix(0, nrow = length(gene_sets), ncol = length(all_genes))
  rownames(binary_matrix) <- paste0("sig.", 1:length(gene_sets))
  colnames(binary_matrix) <- all_genes
  
  for (i in 1:length(gene_sets)) {
    binary_matrix[i, gene_sets[[i]]] <- 1
  }
  
  # Compute cosine similarity
  similarity_matrix <- cosine(t(binary_matrix))
  
  simil.melt <- similarity_matrix %>% reshape2::melt()
  print(
    ggplot(simil.melt, aes(x=Var1,y=Var2,fill=value)) + geom_tile() + theme_bw()+
      scale_fill_gradient2(
        low = "darkblue",
        mid = "white",
        high = "darkred",
        limits = c(0, 1),
        midpoint = 0.5
        
      ) 
  )
  return(gene_sets)
}
