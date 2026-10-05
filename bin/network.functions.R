#! /usr/bin/env Rscript



paper_sources <- c(
  
  "ahn_B"           = 'Ahn et al.',   
  "ahn_A"           = 'Ahn et al.',  
  
  "andersen_class_2"            = 'Andersen et al.',   
  "andersen_class_1"   = 'Andersen et al.',
  
  "song_s100p_pos"              = 'Song et al.',   
  "song_spp1_pos"               = 'Song et al.',   
  
  "martinserrano_Desert-like" = 'Martin-Serrano et al.'        ,
  "martinserrano_Tumor Classical"  = 'Martin-Serrano et al.'  ,
  "martinserrano_Hepati Stem-like"    = 'Martin-Serrano et al.',
  "martinserrano_Inflammatory Stroma"= 'Martin-Serrano et al.',
  "martinserrano_Immune Classical"    = 'Martin-Serrano et al.',
  
  "vijay_bilineage"       = 'Vijay et al.',     
  "vijay_Hepatocyte"       = 'Vijay et al.',     
  "vijay_ductal"          = 'Vijay et al.',     
  "vijay_EMT"             = 'Vijay et al.',     
  
  "orourke_LS"         = 'ORourke et al.',        
  "orourke_RP"                   = 'ORourke et al.',
  
  "montal_mesenchymal"      = 'Montal et al.',   
  "montal_metabolic"             = 'Montal et al.',   
  "montal_proliferation"       = 'Montal et al.',    
  "montal_immune"           = 'Montal et al.',        
  
  "nepal_subtype2"        = 'Nepal et al.',      
  "nepal_subtype1"            = 'Nepal et al.',  
  "nepal_subtype3"           = 'Nepal et al.',  
  "sia_Proliferation"       = 'Sia et al.',    
  "sia_Inflammation"            = 'Sia et al.',   
  "fan_C1"                     = 'Fan et al.',   
  "fan_C2"   = 'Fan et al.'
)

short_name_key <- c(
  "martinserrano_Desert-like" = 'Desert-like'        ,
  "martinserrano_Tumor Classical"  = 'Tumor Classical'  ,
  "martinserrano_Hepati Stem-like"    = 'Hepatic Stem-like',
  "martinserrano_Inflammatory Stroma"= 'Inflammatory Stroma',
  "martinserrano_Immune Classical"   = 'Immune Classical',
  "orourke_LS"         = 'Long Survivor',        
  "orourke_RP"                   = 'Rapid Progressor',
  "vijay_bilineage"       = 'Bilineage',     
  "vijay_Hepatocyte"       = 'Hepatocyte',     
  "vijay_ductal"          = 'Ductal',     
  "vijay_EMT"             = 'EMT',     
  
  
  "montal_mesenchymal"      = 'Mesenchymal',   
  "montal_metabolic"             = 'Metabolic',   
  "montal_proliferation"       = 'Proliferation',    
  "montal_immune"           = 'Immune',        
  "song_s100p_pos"              = 'S100P+',   
  "song_spp1_pos"               = 'SPP1+',   
  "andersen_class_2"            = 'Andersen-2',   
  "andersen_class_1"   = 'Andersen-1',
  "ahn_B"           = 'Ahn-B',   
  "ahn_A"           = 'Ahn-A',  
  "nepal_subtype2"        = 'Nepal-1',      
  "nepal_subtype1"            = 'Nepal-2',  
  "nepal_subtype3"           = 'Nepal-3',  
  "sia_Proliferation"       = 'Sia-Proliferation',    
  "sia_Inflammation"            = 'Sia-Inflammation',   
  "fan_C1"                     = 'Fan-Mesenchymal-C1',   
  "fan_C2"   = 'Fan-Metabolic-C2'
)


assign_consensus_to_calls <- function(labels_df, final.membership){
  #replace class calls with consensus class
  
  labels_df.joined <- labels_df[,-1]
  
  labels_df.joined[] <- lapply(names(labels_df.joined), function(col) {
    paste(col, labels_df.joined[[col]], sep = "_")
  })
  
  labels_df.joined <- as.matrix(labels_df.joined)
  
  cluster_membership.df <- tibble(final.membership)
  cluster_membership.df$signature <- names(final.membership)
  
  lookup_vec <- setNames(cluster_membership.df$final.membership, cluster_membership.df$signature)
  
  mat_num <- matrix(lookup_vec[labels_df.joined], 
                    nrow = nrow(labels_df.joined), 
                    ncol = ncol(labels_df.joined))
  
  
  dt <- as.data.table(mat_num)
  names(dt) <- names(labels_df[,-1])
  rownames(dt) <- labels_df$Sample
  
  dt <- dt %>%
    dplyr::select(where(~ mean(is.na(.)) < 1))
  
  return(dt)
}

assign_classes_to_clusters <- function(final.adj_classed, final.inflation){
  final.g_classed <- graph_from_adjacency_matrix(final.adj_classed, mode = "undirected", weighted = TRUE)
  
  final.mcl_result_classed <- mcl(final.adj_classed, inflation = final.inflation, addLoops = TRUE)
  
  final.membership <- final.mcl_result_classed$Cluster
  names(final.membership) <- rownames(final.adj_classed)
  
  return(final.membership)
}

explore_optimization_space <- function(adj_filtered, inflation_values =  seq(1.5, 20, by = 0.5)){
  
  # Store results
  results <- list()
  g <- graph_from_adjacency_matrix(adj_filtered, mode = "undirected", weighted = TRUE)
  subtypes <- rownames(adj_filtered)
  
  for (inflation in inflation_values) {
    message("\n---------------\nRunning MCL with inflation = ", inflation)
    
    # Run MCL
    mcl_result <- mcl(adj_filtered, inflation = inflation, addLoops = TRUE)
    
    cluster_count <- length(unique(mcl_result$Cluster))
    
    # Extract cluster membership
    membership <- mcl_result$Cluster
    
    # Silhouette score
    # Convert adjacency to distance matrix
    dist_matrix <- as.dist(1 - adj_filtered)  # assuming similarity matrix
    # Compute silhouette
    sil <- silhouette(membership, dist_matrix)
    
    dropped = 0
    
    if ( any(membership <= 0) ) {
      
      message("\nSome signatures have no class.")
      kept_subtypes <- subtypes
      max_iteration = 0
      
      while(any(membership <= 0) ){
        
        message("Dropping ", length(kept_subtypes[membership == 0]), '...')
        dropped = dropped + length(kept_subtypes[membership == 0])
        classed_nodes = kept_subtypes[membership != 0]
        
        # Filter adjacency matrix
        adj_classed <- adj_filtered[classed_nodes, classed_nodes]
        g_classed <- graph_from_adjacency_matrix(adj_classed, mode = "undirected", weighted = TRUE)
        
        # Run MCL
        mcl_result_classed <- mcl(adj_classed, inflation = inflation, addLoops = TRUE)
        
        cluster_count <- length(unique(mcl_result_classed$Cluster))
        
        # Extract cluster membership
        membership <- mcl_result_classed$Cluster
        kept_subtypes <- rownames(adj_classed)
        max_iteration = max_iteration + 1
        
        
        if (max_iteration >= 20) {
          warning("Maximum number of iterations reached. Exiting loop.")
          break
        }
        
        
      }
      
      # Silhouette score
      # Convert adjacency to distance matrix
      dist_matrix_classed <- as.dist(1 - adj_classed)  # assuming similarity matrix
      # Compute silhouette
      sil <- silhouette(membership, dist_matrix_classed)
      
      mod_score <- modularity(g_classed, membership)
      avg_sil <- mean(sil[, 3])
      
      
    } else if (sum(membership) == length(membership)) {
      message("Only one cluster.")
      
      mod_score <- modularity(g, membership)
      avg_sil = NA
      
    } else {
      
      mod_score <- modularity(g, membership)
      avg_sil <- mean(sil[, 3])
      
    }
    
    results <- rbind(results, 
                     data.frame(inflation = inflation, 
                                clusters = cluster_count,
                                modularity = mod_score,
                                silhouette = avg_sil,
                                dropped = dropped
                     ))
  }
  
  
  
  
  # ideal point
  ideal <- c(max(results$modularity, na.rm = T), max(results$silhouette, na.rm = T))
  
  # distance to ideal
  dist <- sqrt((results$modularity - ideal[1])^2 + (results$silhouette - ideal[2])^2)
  
  # row closest to ideal
  final.inflation = results$inflation[which.min(dist) ]
  
  res.melt <- results %>% reshape2::melt(id.vars=c('inflation'))
  
  network_optimization.plot <- 
    ggplot(res.melt, aes(x=inflation, y=value)) + 
    geom_vline(xintercept = final.inflation, color='grey', size=2) +
    geom_line() + geom_point() + facet_grid(variable~., scales='free', switch = "y") +
    theme_bw(base_size = 10)+
    theme(panel.grid = element_blank(), plot.background = element_blank(), strip.background = element_blank(),
          
          panel.background = element_blank(),  axis.title.y = element_blank())  
  
  print(network_optimization.plot)
  output = list(
    'inflation' = final.inflation,
    'inflation.plot' = network_optimization.plot
  )
  
  return(output)
  
}

filter_for_inflation <- function(adj_filtered, final.inflation){
  
  mcl_result <- mcl(adj_filtered, addLoops = TRUE, inflation = final.inflation) 
  interim.membership <- mcl_result$Cluster
  
  subtypes <- rownames(adj_filtered)
  
  final.classed_nodes = subtypes[interim.membership != 0]
  
  # Filter adjacency matrix
  final.adj_classed <- adj_filtered[final.classed_nodes, final.classed_nodes]
  
  dropped.classes <- setdiff(rownames(adj_filtered), rownames(final.adj_classed))
  message('dropping ', length(dropped.classes), ' classes after inflation adjustment: ',
          paste(dropped.classes, collapse = ", "),'\n' )
  
  message(length(final.classed_nodes), ' classes remain: ',
          paste(final.classed_nodes, collapse = ", ") )
  
  return(final.adj_classed)
}

make_adjacency_matrix <- function(labels_df){
  subtype_sets <- list()
  for (col in names(labels_df)[-1]) {
    for (subtype in unique(labels_df[[col]])) {
      subtype_sets[[paste0(col, "_", subtype)]] <- 
        labels_df$Sample[labels_df[[col]] == subtype]
    }
  }
  
  subtypes <- names(subtype_sets)
  n <- length(subtypes)
  jaccard_matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(subtypes, subtypes))
  
  for (i in 1:n) {
    for (j in i:n) {
      a <- subtype_sets[[i]]
      b <- subtype_sets[[j]]
      jaccard <- length(intersect(a, b)) / length(union(a, b))
      jaccard_matrix[i, j] <- jaccard
      jaccard_matrix[j, i] <- jaccard
    }
  }
  
  adj <- jaccard_matrix
  diag(adj) <- 0  
  
  return(adj)
}

make_consensus_cohort <- function(labels_df, dt, treshold = 0.2){
  
  dt_clean <- dt %>%
    dplyr::select(where(~ mean(is.na(.)) <= treshold))
  
  all_1 <- dt_clean[, rowSums(.SD == 1) == length(dt_clean), .SDcols = names(dt_clean)]
  all_2 <- dt_clean[, rowSums(.SD == 2) == length(dt_clean), .SDcols = names(dt_clean)]
  
  dt_clean$tumour.id <- labels_df$Sample
  dt_subset <- dt_clean[all_1 | all_2]
  
  cms_labels <- unlist(c(dt_subset[,1]))
  names(cms_labels) <- dt_subset$tumour.id
  
  message(length(cms_labels), ' samples were consistantly classified across all classes')
  
  CMSO <- setdiff(dt_clean$tumour.id, names(cms_labels))
  
  message(length(CMSO), ' samples were inconsistantly classified')
  
  return(cms_labels)
}

prep_expression_data <- function(basedir = '~/Documents/scripts/github/btc-ccs-data/', threshold = 0.95, subset=NULL){
  
  expr_matrix_raw <- fread(file.path(basedir,"/results/rna.tpm.ff.txt")) 
  expr_matrix <- as.matrix(expr_matrix_raw[,-1])
  rownames(expr_matrix) <- expr_matrix_raw$gene_name
  
  if(is.null(subset)){
    subset = colnames(expr_matrix)
  } else {
    subset = intersect(subset, colnames(expr_matrix))
  }
  
  expr_matrix_raw <- expr_matrix_raw[,..subset]
  
  row_proportions <- rowMeans(expr_matrix > 0)
  filtered_mat <- expr_matrix[row_proportions > threshold, ]
  
  gene_set <- fread(file.path(basedir,'/source_data/fgs.txt'))
  
  gene_list <- gene_set %>% 
    filter( ( biotype=="protein_coding" & Length >= 500 & Chr %in% c(1:23) ) ) %>% #
    dplyr::select(gene_name)
  
  these.genes <- intersect(gene_list$gene_name, rownames(filtered_mat))
  filtered_mat <- filtered_mat[these.genes,]
  
  return(filtered_mat)
}

remove_disconnected_classes <- function(labels_df, adj, cutoff = 0.02){
  
  keep_nodes <- which(rowSums(adj) >= nrow(labels_df) * cutoff)
  adj_filtered <- adj[keep_nodes, keep_nodes]
  
  message('dropping ', length(setdiff(rownames(adj), rownames(adj_filtered))), ' class for low connectivity: ',
          paste(setdiff(rownames(adj), rownames(adj_filtered)), collapse = ", "))
  return(adj_filtered)
  
}

remove_unassigned_classifiers <- function(labels_df.raw, classifieds.n, cutoff=0.5){
  #filter.predictions
  
  labels_df <- labels_df.raw %>% dplyr::select(starts_with("Sample") | ends_with('.prediction'))
  names(labels_df) <- gsub('.prediction','',names(labels_df) )
  
  low.assignment.classifiers <- classifieds.n %>% filter(is.na(class) & percent > cutoff) %>% 
    ungroup() %>% dplyr::select(variable) %>% c() %>% unlist()
  
  message('removing ',length(low.assignment.classifiers), ' classifiers for low assignment rates: ',
          paste(low.assignment.classifiers, collapse = ", "),'\n' 
  )
  
  labels_df <- labels_df  %>% dplyr::select( -any_of(low.assignment.classifiers ) )
  
  message(ncol(labels_df) -1, ' classifiers remain.')
  
  return(labels_df)
}

tally_class_assignments <- function(labels_df.raw){
  # how many to remove
  
  classifieds <- labels_df.raw %>% dplyr::select(starts_with("Sample") | ends_with('.class'))
  
  message('processing ',
          nrow(classifieds),
          ' samples across ',
          length(unique(unlist(classifieds[,-1]))) ,
          ' classes from ',
          ncol(classifieds) - 1,
          ' classifiers.' 
  )
  
  classifieds.melt <- classifieds %>% reshape2::melt(id.vars=c('Sample')) 
  classifieds.melt$variable <- gsub(".class","", classifieds.melt$variable)
  classifieds.melt$class <- paste0(classifieds.melt$variable, "_", classifieds.melt$value)
  classifieds.melt$class[is.na(classifieds.melt$value)] <- NA
  
  classifieds.n <- classifieds.melt  %>% group_by(variable, value, class) %>% tally()
  classifieds.n$percent = classifieds.n$n / nrow(classifieds)
  return(classifieds.n)
}

train_tsp_on_consensus_cohort <- function(filtered_mat, cms_labels, data.dir, save.train.test = F){
  
  labeled_samples <- intersect(colnames(filtered_mat), names(cms_labels))
  X.raw <- filtered_mat[, labeled_samples]
  
  ## this currently does nothing but could be adapted to subset to variable genes
  gene_variances <- apply(X.raw, 1, var)
  top_genes <- names(sort(gene_variances, decreasing = TRUE))
  expr_filtered <- X.raw[top_genes, ]
  
  X <- t(expr_filtered)
  
  y <- cms_labels[labeled_samples]
  
  message('full set has ', length(y), ' samples')
  
  train_idx <- sample(seq_along(y), 4 * length(y) / 5)
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  
  message('training set has ', length(y_train), ' samples')
  
  X_valid <- X[-train_idx, ]
  y_valid <- y[-train_idx]
  
  if(save.train.test){
    tsp.training.samples <- list('training_set'=y_train, 'test_set'=y_valid)
    saveRDS(tsp.training.samples, paste0(data.dir,'results/tsp.training.samples.rds'))
  }
  
  #### TSP model ####
  
  expr_data <- t(X_train)
  
  trainingGroupNum <- as.numeric(y_train) - 1
  
  table(trainingGroupNum)
  
  classifier <- SWAP.KTSP.Train(inputMat = expr_data, 
                                phenoGroup = factor(trainingGroupNum), 
                                krange = c(10:16))
  
  trainingPrediction <- SWAP.KTSP.Classify(t(X_valid), classifier)
  
  acc.table <- table(Predicted = trainingPrediction, Actual = y_valid)
  print(acc.table)
  correct.n = acc.table[1,1] + acc.table[2,2]
  
  print(BinomCI(x = correct.n, n = sum(acc.table), conf.level = 0.95, 
                sides = c("two.sided"), method = c("wilson")))
  
  
  #https://www.genenames.org/download/statistics-and-files/
  HGNC <- read.delim(file.path(data.dir,"/source_data/hgnc_complete_set.txt"), sep = "\t", stringsAsFactors = FALSE)
  
  classifier[['aliases']] <- HGNC %>% filter(symbol %in% c(unlist(data.frame(classifier$TSPs)))) %>% dplyr::select(symbol, alias_symbol, prev_symbol)
  
  return(classifier)
}

