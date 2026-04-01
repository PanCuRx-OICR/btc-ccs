
make_cn_tpm_table <- function(this_gene, BTC.cnv.matrix, cnv.gene_list, all_rna, rna.gene_list){
  cnv_named <-  cbind.data.frame(
    'tumour' = names(BTC.cnv.matrix),
    'cnv' = unlist(BTC.cnv.matrix[cnv.gene_list == this_gene,]),
    'cnv.rank'= rank(BTC.cnv.matrix[cnv.gene_list == this_gene,])
  )
  
  tpm_named <-  cbind.data.frame(
    'tumour' = names(all_rna),
    'tpm'= unlist(all_rna[rna.gene_list == this_gene,]),
    'tpm.rank'= rank(all_rna[rna.gene_list == this_gene,])
  )
  
  cnv.tpm <- inner_join(cnv_named, tpm_named, by=c('tumour'='tumour'))
  return(cnv.tpm)
}


get_cn_tpm_correlations <- function(this_gene, BTC.cnv.matrix, cnv.gene_list, all_rna, rna.gene_list ){
  
  cnv.tpm <- make_cn_tpm_table(this_gene, BTC.cnv.matrix, cnv.gene_list, all_rna, rna.gene_list)
  
  this.lm <- lm( tpm ~  cnv , cnv.tpm)
  this.lm <- summary(this.lm)
  
  if(Inf %in% log(cnv.tpm$cnv) | -Inf %in% log(cnv.tpm$cnv)){
    log.lm.p.coefficient = NA
    log.lm.r = NA
  } else {
    this.log.lm <- lm( tpm ~  log(cnv) , cnv.tpm)
    this.log.lm <- summary(this.log.lm)
    log.lm.p.coefficient = this.log.lm$coefficients[[8]]
    log.lm.r = this.log.lm$adj.r.squared
  }
  
  this.corr <- cor.test(x=cnv.tpm$cnv, y=cnv.tpm$tpm, method = 'spearman')
  
  this.tmp.rel <- cbind.data.frame('gene'=this_gene, 
                                   'lm.p'=this.lm$coefficients[[8]], 
                                   'lm.r2'=this.lm$adj.r.squared, 
                                   'log.lm.p'= log.lm.p.coefficient,
                                   'log.lm.r2'= log.lm.r,
                                   'sr.p' = this.corr$p.value,
                                   'sr.rho' = this.corr$estimate
  )
  return(this.tmp.rel)
}


filter_and_annotate_centromeres <- function(scores.gistic, centromeres, manual_filter_file, MIN_PEAK_AMPLITUDE = 10, MIN_SIZE_FROM_CENTROMERE = 5000000){
  REMOVE_THESE_CHROMOSOMES = c("X","Y")
  
  manually_remove_gistic <- fread(manual_filter_file)
  
  for(this_row_index in c(1:nrow(scores.gistic))){
    this_chr <- scores.gistic$Chromosome[this_row_index]
    
    if(scores.gistic$minus_log_qvalue[this_row_index] > MIN_PEAK_AMPLITUDE ){
      
      
      this_centromere_pos <- centromeres$V2[centromeres$V1 == this_chr]
      
      cat("checking peak on ", this_chr, "\n")
      
      if(abs(scores.gistic$Start[this_row_index] - this_centromere_pos) < MIN_SIZE_FROM_CENTROMERE){
        scores.gistic$minus_log_qvalue[this_row_index] <- 0
        cat("removing peak on ", this_chr, " within ", abs(scores.gistic$Start[this_row_index] - this_centromere_pos), "\n")
        
      }
      
     
      
    }
    for(manual_region_index in c(1:nrow(manually_remove_gistic))){
      if(manually_remove_gistic$V2[manual_region_index] == paste0('chr',this_chr)){
        if(scores.gistic$Start[this_row_index] <= (manually_remove_gistic$V4[manual_region_index] ) & 
           scores.gistic$End[this_row_index] >= (manually_remove_gistic$V3[manual_region_index] )
        ){
          scores.gistic$minus_log_qvalue[this_row_index] <- 0
        }
      }
    }
    
    
  }
  
  centromeres <- centromeres[centromeres$V1 %ni% REMOVE_THESE_CHROMOSOMES,]
  centromeres$holder <- NA
  centromeres$Type <- "Cent"
  centromeres_rearr <- centromeres[,c("Type", "V1", "V2", "V2", "holder", "holder", "holder", "holder", "holder")]
  names(centromeres_rearr) <- names(scores.gistic)
  
  scores.gistic <- rbind.data.frame(scores.gistic, centromeres_rearr)
  scores.gistic$cent[scores.gistic$Type == "Cent"] <- T
  return(scores.gistic)
}

collate_barebones_vcf_write <- function(sample.list, inpath, outpath='~/Documents', variant_types = c("indel","snv","dbs")){
  first=T
  for(tumor in sample.list$V2){
    
    cat(tumor,"\n")
    for(variant_type in variant_types){ 
      file_name <- paste0(inpath, '/',tumor, '/pyclone/',tumor,'.',variant_type,'.barebones.vcf')
      
      if(file.exists(file_name)){
        
        these_mutations = fread(file_name)
        these_mutations$patient <- tumor
        
        file_name_info <- paste0(inpath, '/',tumor, '/pyclone/',tumor,'.',variant_type,'.info.txt')
        these_mutations_info = fread(file_name_info)
        these_mutations_info$patient <- tumor
        
        
        if(first){
          all_mutations <- these_mutations 
          all_info <- these_mutations_info 
          first=F}
        else{
          all_mutations <- rbind.data.frame(all_mutations, these_mutations)
          all_info <- rbind.data.frame(all_info, these_mutations_info)
          
        }
      } else {
        cat(file_name,' not found\n')
      }
    }
  }
  
  write.table(
    all_mutations,
    file = paste0(outpath,"/all_mutations.txt"),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  
  write.table(
    all_info,
    file = paste0(outpath,"/all_info.txt"),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
}

get_maxmin_from_splits <- function(cn_string_vec, maxmin) {
  max_vec <- c()
  #need to convert to character otherwise vecs where all are numeric come in as numeric
  #and strsplit won't split numerics
  cn_string_vec <- as.character(cn_string_vec)
  for(this_cn_string_index in c(1:length(cn_string_vec))){
    
    this_num_vec <- as.numeric(unlist(strsplit(as.character(cn_string_vec[this_cn_string_index]), "\\|")))
    if(maxmin=='AMP'){
      this_max <- max(this_num_vec)
    } else {
      ## == 'DEL'
      this_max <- min(this_num_vec)
      
    }
    max_vec <- c(max_vec, this_max)
  }
  return(max_vec)
}

get_dnds_sig_genes <- function(sel_cv, dnds_cutoff = 0.1, INDELS_BOOLEAN = T){
  
  if(INDELS_BOOLEAN){
    #only with indels
    signif_genes = sel_cv[sel_cv$qglobal_cv<dnds_cutoff, c("gene_name","qglobal_cv")]
    rownames(signif_genes) = NULL
  }else{
    
    signif_genes = sel_cv[sel_cv$qallsubs_cv<dnds_cutoff, c("gene_name","qallsubs_cv")]
    rownames(signif_genes) = NULL
  }
  return(signif_genes)
}


pull_clonality <- function(these_dnds_annotation, work_dir_from_local='/Volumes/pcsi/users/fbeaudry/bunch_o_vcfs/'){
  require(tidyr)
  first=T
  
  for(tumour_sample_id in unique(these_dnds_annotation$sampleID)){
    cat(tumour_sample_id, "\n")
    warning_printed = F
    mutation_timer_file_name <- paste0(work_dir_from_local,tumour_sample_id,'.mutationTimeR.annotated.txt')
    pyclone_file_name <- paste0(work_dir_from_local,tumour_sample_id,'.pyclone.tsv')
    
    if(file.exists(mutation_timer_file_name) & file.exists(pyclone_file_name)){
      
      mutationTimeR.annotated <- fread(mutation_timer_file_name)
      names(mutationTimeR.annotated)[1] <- "CHROM"
      
      pyclone <- fread(pyclone_file_name)
      pyclone <- separate(data = pyclone, col = mutation_id, into = c("chr", "pos", "sub"), sep = ':')
      pyclone$pos <- round(as.numeric(pyclone$pos))
      #unique(pyclone[,c("cluster_id", "cellular_prevalence")])
      
      alterations_in_this_sample <- these_dnds_annotation[these_dnds_annotation$sampleID == tumour_sample_id,]
      
      for(this_alteration_index in c(1:nrow(alterations_in_this_sample))){
        #cat(this_alteration_index, " ")
        #current data only includes SNVs
        if(alterations_in_this_sample$impact[this_alteration_index] != 'no-SNV'){
          
          this_chr <- paste0('chr',alterations_in_this_sample$chr[this_alteration_index])
          this_pos <- alterations_in_this_sample$pos[this_alteration_index]
          this_clonality <- mutationTimeR.annotated[mutationTimeR.annotated$POS == this_pos & mutationTimeR.annotated$CHROM == this_chr,]
          this_subclonality <- pyclone[pyclone$pos == this_pos & pyclone$chr == this_chr,]
          if(nrow(this_subclonality) > 0){
            cluster_id = this_subclonality$cluster_id
            cellular_prevalence = this_subclonality$cellular_prevalence
          } else {
            cluster_id = NA
            cellular_prevalence = NA
          }
          
          clonality_table_tmp <- cbind.data.frame('tumour_sample_id' = tumour_sample_id, 
                                                  'chr' = this_chr,
                                                  'pos' = this_pos,
                                                  'gene'=alterations_in_this_sample$gene[this_alteration_index], 
                                                  'aachange' = alterations_in_this_sample$aachange[this_alteration_index],
                                                  'impact' = alterations_in_this_sample$impact[this_alteration_index],
                                                  'clonality' = this_clonality$CLS,
                                                  'cluster_id' = cluster_id,
                                                  'cellular_prevalence' = cellular_prevalence
          )
          
          if(first){clonality_table <- clonality_table_tmp; first=F}else{clonality_table <- rbind.data.frame(clonality_table, clonality_table_tmp)}
        } else {
          if(warning_printed == F){
            cat('skipping non-SNV variants for ',tumour_sample_id,'\n')
            warning_printed = T
          }
        }
      }
    } else {
      cat(mutation_timer_file_name, " or ", pyclone_file_name, " not found\n")
    }
  }
  return(clonality_table)
}

split_at_8th <- function(x) {
  # Extract the first part (up to the 8th character)
  part1 <- substr(x, 1, 8)
  return(unlist(part1))
}
