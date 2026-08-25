#! /usr/bin/env Rscript

assemble_duct_markers <- function(data.dir, gene_set){
  sampaziotis_2021_tableS1 <- fread(file.path(data.dir,'/source_data/sampaziotis_2021_tableS1.txt'), header=T)
  sampaziotis_2021_tableS1 <- left_join(sampaziotis_2021_tableS1, gene_set, by=c('Gene_name'='gene_name'), relationship =  "many-to-many")
  sampaziotis_gene_sets <- sampaziotis_2021_tableS1 %>% dplyr::select(Region, ensid) %>% filter(!is.na(ensid))
  sampaziotis_gene_sets$Region <- paste0("sampaziotis.",sampaziotis_gene_sets$Region)
  
  rimland_2021_tableS2 <- fread(file.path(data.dir,'/source_data/rimland_2021_tableS2.txt'), header=T)
  rimland_gene_sets <- rimland_2021_tableS2 %>% dplyr::select(V3, `ENSG ID`)
  rimland_gene_sets$V3 <- paste0("rimland.",rimland_gene_sets$V3)
  names(rimland_gene_sets) <- names(sampaziotis_gene_sets)
  
  andrews2024_tableS10_sn <- fread(file.path(data.dir,'/source_data/andrews_tableS10_sn.celltypemarkers.txt'), header=T)
  andrews2024_tableS10_sn <- andrews2024_tableS10_sn %>% filter(avg_logFC > 0)
  andrews2024_tableS10_sn <- left_join(andrews2024_tableS10_sn, gene_set, by=c('gene'='gene_name'), relationship =  "many-to-many")
  andrews2024sn_gene_sets <- andrews2024_tableS10_sn %>% dplyr::select(cluster, ensid) %>% filter(!is.na(ensid)) %>%
    filter(!apply(., 1, function(row) any(str_detect(row, "Doublet"))))
  
  names(andrews2024sn_gene_sets) <- names(sampaziotis_gene_sets)
  
  duct.gene.set <- rbind.data.frame(sampaziotis_gene_sets, rimland_gene_sets, andrews2024sn_gene_sets)
  return(duct.gene.set)
}

full_correlation_matrix <- function(H_matrix_ordered_norm, order_vec, component_n){
  
  cor_with_p <- function(x) {
    n <- ncol(x)
    cor_mat <- matrix(NA, n, n)
    p_mat <- matrix(NA, n, n)
    colnames(cor_mat) <- colnames(x)
    rownames(cor_mat) <- colnames(x)
    colnames(p_mat) <- colnames(x)
    rownames(p_mat) <- colnames(x)
    
    for (i in 1:n) {
      for (j in 1:n) {
        test <- cor.test(x[, i], x[, j])
        cor_mat[i, j] <- test$estimate
        p_mat[i, j] <- test$p.value
      }
    }
    
    cor_mat[p_mat > 0.05] <- NA  
    return(cor_mat)
  }
  
  cor_matrix <- cor_with_p(t(H_matrix_ordered_norm))
  
  cor_matrix <- cor_matrix[order_vec, order_vec]
  
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] <- NA
  
  cor_df <- data.frame(cor_matrix)
  names(cor_df) <- gsub('X','sig',names(cor_df))
  cor_df$sig.id <- paste0('sig',rownames(cor_matrix))
  
  cor_melt <- reshape2::melt(cor_df, id.vars='sig.id')
  cor_melt_filt <- cor_melt %>% filter(!is.na(value))
  cor_melt_filt <- left_join(data.frame('sig.id'=paste0('sig',c(1:component_n))), cor_melt_filt, by='sig.id') 
  cor_melt_filt$cor.bin <- NA
  cor_melt_filt$cor.bin[cor_melt_filt$value > 0 ] <- 'up'
  cor_melt_filt$cor.bin[cor_melt_filt$value < 0 ] <- 'down'
  return(cor_melt_filt)
}

get_ccs_enrichment <- function(W_matrix, data.dir, component_n, y_order.convert){
  tsp.genes <- fread(file.path(data.dir,'/results/tsp.genes.txt'))
  names(tsp.genes) <- gsub('M','C',names(tsp.genes) )
  tsp.genes$gene <- 'gene'
  tsp.genes <- melt(tsp.genes, id.vars = 'gene') %>% 
    left_join(gene_set, by=c('value'='gene_name')) %>% 
    dplyr::select(variable, ensid)
  
  first.loop = T
  cms.gsea <- NULL
  for(this.sig in c(1:ncol(W_matrix)) ){
    
    cat('--\n\n',this.sig,"\n")
    this.sig.w <- W_matrix[,this.sig]
    names(this.sig.w) <- rownames(W_matrix)
    these_ranked_genes <- sort(this.sig.w, decreasing = TRUE)
    
    gsea_result <- GSEA(geneList = these_ranked_genes,
                        TERM2GENE = tsp.genes,
                        verbose = T,
                        
                        scoreType = "pos", minGSSize = 10,
                        pvalueCutoff = 1)
    
    
    gseaplot2(gsea_result, geneSetID = "CCS-A")
    gseaplot2(gsea_result, geneSetID = "CCS-B")
    
    if(nrow(gsea_result@result) > 0){
      gsea.df <- gsea_result@result[,c(1:10)]
      gsea.df$sig.id <- paste0('sig',this.sig)
      if(first.loop){cms.gsea <- gsea.df; first.loop=F}else{cms.gsea <- rbind.data.frame(cms.gsea,gsea.df)}
      
    }
  }
  
  if(!is.null(cms.gsea)){
    cms.gsea <- left_join(data.frame('sig.id'=paste0('sig',c(1:component_n))), cms.gsea, by=c('sig.id'='sig.id')) %>% 
      left_join( y_order.convert, by=c('sig.id'='og'))
    
  }
  cms.gsea <- cms.gsea %>% arrange(desc(p.adjust))
  
  return(cms.gsea)
}

get_cibersort_correlation <- function(cibersort_path, H_matrix_ordered_norm.melt, component_n){
  
  LBR.cibersort <- fread(cibersort_path) %>% melt(id.vars='cell_type')
  
  first.loop = T
  ## loop across components/signatures
  for(this.sig in unique(H_matrix_ordered_norm.melt$sig.id)){
    
    ## filter for weight higher than median (zero, after scale)
    this.sig.w <- H_matrix_ordered_norm.melt %>% 
      filter(sig.id == this.sig  ) %>% 
      left_join(LBR.cibersort, by=c('variable'='variable')) %>%
      filter(!is.na(cell_type))
    
    for(this.cell.type in unique(this.sig.w$cell_type)){
      
      this.cell.type.w <- this.sig.w %>% 
        filter(cell_type == this.cell.type  ) 
      
      ## calculate correaltion between sig weight and cellularity
      this.cor <- cor.test(this.cell.type.w$value.x, this.cell.type.w$value.y)
      
      ## make row
      this.cor <- cbind.data.frame(
        'sig.id'=this.sig,
        'cell.type'=this.cell.type,
        'cor'=this.cor$estimate,
        'cor.p'=this.cor$p.value
      )
      
      if(first.loop){immune.cor <- this.cor; first.loop=F}else{immune.cor <- rbind.data.frame(immune.cor,this.cor)}
    }
  }
  
  
  immune.cor.filt <- left_join(data.frame('sig.id'=paste0('sig',c(1:component_n))), 
                               immune.cor, by=c('sig.id'='sig.id')) 
  
  return(immune.cor.filt)
}

get_clustering_order <- function(cellularity.sigs.mat){
  ph <- pheatmap(
    t(cellularity.sigs.mat),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    clustering_distance_cols = "euclidean",
    clustering_method = "complete", 
    color = colorRampPalette(c("darkblue","white", "darkgreen"))(100),
    
    show_colnames = TRUE,
    fontsize = 8,
    border_color = "white"
  )
  
  x_order <- ph$tree_col$labels[ph$tree_col$order]
  
  heat.cluster <- data.frame(
    sample = x_order,
    order  = seq_along(x_order),
    stringsAsFactors = FALSE
  )
  return(heat.cluster)
}

get_cophenetic_coef <- function(cellularity.sigs.mat){
  hc_rows <- hclust(dist(cellularity.sigs.mat))
  
  # Calculate the cophenetic correlation coefficient
  correlation <- cor(dist(cellularity.sigs.mat), cophenetic(hc_rows))
  return(correlation)
}

get_mscd <- function(filtered_genes_cts, nmf_results){
  ## Reconstruct and compute cosine distance (MSCD)
  V_rec <- fitted(nmf_results@fit)
  V_rec_norm <- t(apply(V_rec, 1, function(x) x / sum(x)))
  
  og_norm <- t(apply(filtered_genes_cts, 1, function(x) x / sum(x)))
  sample_cosine_distances <- proxy::dist(og_norm, V_rec_norm, method = "cosine")
  return(mean(as.vector(sample_cosine_distances)))
}

get_signature_stats <- function(this.sig, W_matrix, H_matrix_ordered_norm.melt, verbose=F){
  ## filter for weight higher than median (zero, after scale)
  this.sig.w <- H_matrix_ordered_norm.melt %>% 
    filter(sig.id == this.sig & !is.na(cellularity) & value > 0 ) %>%
    dplyr::select(value, cellularity)
  
  ## calculate correaltion between sig weight and cellularity
  this.cor <- cor.test(this.sig.w$value, this.sig.w$cellularity)
  
  ## annotate sig cell type
  this.sig.col <- paste0('V', gsub('sig','',this.sig))
  
  this.W_matrix <- W_matrix %>% dplyr::select(all_of(this.sig.col)) 
  names(this.W_matrix) <- 'w'
  
  ## pull 100 top genes in component for annotation
  these.genes <- row.names(this.W_matrix %>% arrange(desc(w)) %>% slice_head(n=100))
  
  sig.profile <- gost(query = these.genes, organism = "hsapiens", sources = 'HPA', significant = FALSE)
  
  ## if successfully matched to tissue, annotate in table
  if(length(sig.profile$result) > 0){
    
    sig.profile.res <- sig.profile$result
    sig.profile.res$term_name <- gsub('-', ';', sig.profile.res$term_name)
    top.tissue = paste0(unlist(sig.profile.res %>%  slice_head(n = 3) %>% dplyr::select(term_name)),  collapse = "-")
    top.tissue.p = paste0(unlist(sig.profile.res %>%  slice_head(n = 3) %>% dplyr::select(p_value)),  collapse = ";")
    top.tissue.precision = paste0(unlist(sig.profile.res %>%  slice_head(n = 3) %>% dplyr::select(precision)),  collapse = ";")
    top.tissue.recall = paste0(unlist(sig.profile.res %>%  slice_head(n = 3) %>% dplyr::select(recall)),  collapse = ";")
    
    if(verbose){message(this.sig,': ' ,top.tissue,'\n')}
    
    
  } else {
    top.tissue = 'none'
    top.tissue.p=NA
    top.tissue.precision=NA
    top.tissue.recall=NA
  }
  
  ## make row
  this.cor <- cbind.data.frame(
    
    'n'= nrow(this.sig.w),
    'sig.id'=this.sig,
    'cor'=this.cor$estimate,
    'cor.p'=this.cor$p.value,
    'top.tissue'=top.tissue,
    'top.tissue.p'=top.tissue.p,
    'top.tissue.precision'=top.tissue.precision,
    'top.tissue.recall'=top.tissue.recall
  )
  return(this.cor)
}

get_tissue_site <- function(site.info, H_matrix_ordered_norm.melt, compare_by_t = F,  compare_by_f = T){
  
  H.by_site <- left_join(H_matrix_ordered_norm.melt, site.info, by=c('variable'='Sample')) %>%
    filter(!is.na(Site_of_tissue))
  
  site.test = NULL
  for(this.sig in unique(H.by_site$sig.id)){
    
    for(this.site in unique(H.by_site$Site_of_tissue)){
      
      this.site.test <- H.by_site %>% filter(sig.id == this.sig)
      
      if(compare_by_t){
        t.res <- t.test(this.site.test$value[this.site.test$Site_of_tissue == this.site],
                        this.site.test$value[this.site.test$Site_of_tissue != this.site])
        
        this.res <- cbind.data.frame(
          'sig.id'=this.sig,
          'site'=this.site,
          'effect'=abs(t.res$estimate[1]) + abs(t.res$estimate[2]) ,
          'p'=t.res$p.value
        )
      }
      
      if(compare_by_f){
        this.site.test$hit <- NA
        this.site.test$hit[this.site.test$value >= 0] <- 'signature'
        this.site.test$hit[this.site.test$value < 0] <- 'no signature'
        
        this.site.test$site[this.site.test$Site_of_tissue == this.site] <- 'this site'
        this.site.test$site[this.site.test$Site_of_tissue != this.site] <- 'other site'
        
        f.res <- this.site.test %>% dplyr::select(hit, site) %>% table() %>% fisher.test()
        
        this.res <- cbind.data.frame(
          'sig.id' = this.sig,
          'site' = this.site,
          'effect'= f.res$estimate ,
          'p' = f.res$p.value
        )
      }
      
      if(is.null(site.test)){site.test <- this.res} else {site.test <- rbind.data.frame(site.test,this.res)}
      
    }
    
  }
  site.test$effect[is.infinite(site.test$effect)] <- max(site.test$effect[!is.infinite(site.test$effect)])
  site.test$site <- factor(site.test$site, levels=levels(site.info$Site_of_tissue))
  
  return(site.test)
}

pull_tissue_match <- function(nmf.stats, component_n){
  
  tissues.matches <- nmf.stats %>% 
    dplyr::select(sig.id, top.tissue, top.tissue.p, top.tissue.precision, top.tissue.recall )
  
  tissues.matches.melt <- 
    separate(tissues.matches, top.tissue, into = c('top.tissue.1','top.tissue.2','top.tissue.3'), sep = "-") %>% 
    separate(top.tissue.p, into = c('top.tissue.1.p','top.tissue.2.p','top.tissue.3.p'), sep = ";") %>% 
    separate(top.tissue.precision, into = c('top.tissue.1.precision','top.tissue.2.precision','top.tissue.3.precision'), sep = ";") %>% 
    separate(top.tissue.recall, into = c('top.tissue.1.recall','top.tissue.2.recall','top.tissue.3.recall'), sep = ";") %>%
    
    reshape2::melt(id.vars=c('sig.id'))  %>% 
    separate( variable, into = c('top','tissue','hit','var'), sep = "\\.")
  
  tissues.matches.melt$var[is.na(tissues.matches.melt$var)] <- 'tissue'
  
  tissues.matches.df <- reshape2::dcast(data = tissues.matches.melt, formula = sig.id + hit ~ var) 
  
  tissues.matches.df <- separate(tissues.matches.df, tissue, into = c('tissuetype','celltype'), sep = "; ", remove = F)
  tissues.matches.df$celltype <- gsub('\\[\\≥Medium\\]','',tissues.matches.df$celltype)
  tissues.matches.df$celltype <- gsub('\\[\\≥Low\\]','',tissues.matches.df$celltype)
  tissues.matches.df$celltype <- gsub('\\[High]','',tissues.matches.df$celltype)
  
  tissues.matches.df <- left_join(data.frame('sig.id'=paste0('sig',c(1:component_n))), tissues.matches.df, by=c('sig.id'='sig.id')) 
  
  ## some extra cleaning, maybe unnecessary
  tissues.matches.df$tissuetype[tissues.matches.df$tissuetype %in% c('Skin 1','Skin 2')] <- 'Skin'
  tissues.matches.df$celltype <- gsub(" \\(cell body\\)","",tissues.matches.df$celltype)
  
  return(tissues.matches.df)
}

reformat_h_matrix <- function(H_matrix, cohort.info=NULL){
  
  H_matrix_ordered_norm <-   data.frame(t(apply(H_matrix, 1, scale)))
  names(H_matrix_ordered_norm) <- names(H_matrix)
  
  H_matrix_ordered_norm.df <- data.table(H_matrix_ordered_norm)
  H_matrix_ordered_norm.df$sig.id <- paste0('sig', row.names(H_matrix_ordered_norm))
  
  H_matrix_ordered_norm.melt <- melt(H_matrix_ordered_norm.df, id.vars = 'sig.id')
  
  if(!is.null(cohort.info)){
    H_matrix_ordered_norm.melt <- H_matrix_ordered_norm.melt %>% 
      left_join(cohort.info, by=c('variable'='tumour'))
  }
  
  return(H_matrix_ordered_norm.melt)
}

reformat_w_matrix <- function(W_matrix, gene_kb){
  tpm_data <- exp(W_matrix) - 1
  tpm_data <- tpm_data / gene_kb[rownames(W_matrix)]
  scaling_factors <- colSums(tpm_data)
  tpm_data <- (tpm_data/ scaling_factors) * 1e6
  tpm_data <- log(tpm_data + 1)
  return(tpm_data)
}

test_duct_set <- function(W_matrix, duct.gene.set, component_n){
  first.loop = T
  for(this.sig in c(1:ncol(W_matrix)) ){
    
    message('--\n',this.sig)
    this.sig.w <- W_matrix[,this.sig]
    names(this.sig.w) <- rownames(W_matrix)
    these_ranked_genes <- sort(this.sig.w, decreasing = TRUE)
    
    gsea_result <- GSEA(geneList = these_ranked_genes,
                        TERM2GENE = duct.gene.set,
                        scoreType = "pos",
                        pvalueCutoff = 1)
    
    if(nrow(gsea_result@result) > 0){
      gsea.df <- gsea_result@result[,c(1:10)]
      gsea.df$sig.id <- paste0('sig',this.sig)
      if(first.loop){chol.gsea <- gsea.df; first.loop=F}else{chol.gsea <- rbind.data.frame(chol.gsea,gsea.df)}
      
    }
  }
  
  chol.gsea$Study <- 'Andrews et al.' 
  chol.gsea$Study[chol.gsea$Description %in% c('rimland.CBD','rimland.GBD','rimland.PancD')] <- 'Rimland et al.' 
  chol.gsea$Study[chol.gsea$Description %in% c('sampaziotis.CBD','sampaziotis.GB','sampaziotis.IHD')] <- 'Sampaziotis et al.' 
  
  chol.gsea$p.adjust.multi <- chol.gsea$p.adjust * component_n * length(unique(chol.gsea$Description))
  chol.gsea$p.adjust.multi[chol.gsea$p.adjust.multi > 0.05] <- 1
  
  return(chol.gsea)
}

save.nmf.matrices <- function(component_n, gene_n, cluster.dir, data.dir, cohort=''){
  
  H_matrix <- read.csv(file.path(cluster.dir, paste0("h_mat.k",component_n,".",gene_n,"g.csv")), row.names = 1)
  
  H_matrix.print <-   data.frame(apply(H_matrix, 1, scale))
  H_matrix.print <- cbind.data.frame('tumour.id'=names(H_matrix), H_matrix.print)
  names(H_matrix.print) <- gsub('X','sig', names(H_matrix.print))
  
  W_matrix <- read.csv(file.path(cluster.dir,paste0("w_mat.k",component_n,".",gene_n,"g.csv")), row.names= 1)
  
  W_matrix.print <- W_matrix
  W_matrix.print$gene_name <- row.names(W_matrix)
  W_matrix.print <- left_join(W_matrix.print, gene_set, by=c('gene_name'='gene_name'))
  
  write.table( H_matrix,  
               file = file.path(data.dir,'results', paste0(cohort,'.H_matrix.txt')), 
               row.names = T,col.names = T, quote = FALSE, sep = ",", append =F)
  write.table( H_matrix.print,  
               file = file.path(data.dir,'results', paste0(cohort,'.H_matrix.print.txt')), 
               row.names = F,col.names = T, quote = FALSE, sep = "\t", append =F)
  
  write.table( W_matrix,  
               file = file.path(data.dir,'results', paste0(cohort,'.W_matrix.txt')), 
               row.names = T,col.names = T, quote = FALSE, sep = ",", append =F)
  write.table( W_matrix.print,  
               file = file.path(data.dir,'results', paste0(cohort,'.W_matrix.print.txt')), 
               row.names = F,col.names = T, quote = FALSE, sep = "\t", append =F)
}

sort_celltypes <- function(chol.gsea.hits){
  
  unique_sc_celltypes <- unique(chol.gsea.hits$cell[order(chol.gsea.hits$new, decreasing = T)])
  sorted_celltypes <- chol.gsea.hits$cell[order(chol.gsea.hits$new, decreasing = T)]
  
  first_indices <- match(unique_sc_celltypes, sorted_celltypes)
  chol.gsea.hits$cell.f <- factor(chol.gsea.hits$cell, levels = sorted_celltypes[first_indices])
  return(chol.gsea.hits)
}
