'%ni%' <- function(x,y)!('%in%'(x,y))

concatenate_ids <- function(this_req){
  # concatenate IDs
  first_id=T
  for(this_id in c(1:length(this_req$ids))){
    if(first_id){id_list <- this_req$ids[[this_id]]; first_id=F}
    else{id_list <- c(id_list,this_req$ids[[this_id]] )}
    
  }
  id_list <- sort(id_list)
  pasted_id <- paste(id_list, collapse=", ")
  
  return(pasted_id)
}

get_gene_frequencies <- function(variants_by_sample){
  # makes sure donors can only have a mutation once (no two-knockouts, no two samples)
  by_donor <- variants_by_sample %>% group_by(donor, gene) %>% tally()
  by_donor$n <- 1
  variant_frequency <- by_donor %>% group_by(gene) %>% 
    reframe(frequency = sum(n) / length(unique(variants_by_sample$donor)))
  
  variant_frequency$gene_factor    <- factor(variant_frequency$gene,  levels = levels(variants_by_sample$gene_factor))
  return(variant_frequency)
}


get_gene_frequencies_oncokb <- function(variants_by_sample){
  # makes sure donors can only have a mutation once (no two-knockouts, no two samples)
  by_donor <- variants_by_sample %>% group_by(donor, gene, top_oncokb_level) %>% tally()
  by_donor$n <- 1
  variant_frequency <- by_donor %>% group_by(gene, top_oncokb_level) %>% 
    reframe(tally=sum(n), frequency = sum(n) / length(unique(variants_by_sample$donor)))
  
  variant_frequency$gene_factor    <- factor(variant_frequency$gene,  levels = levels(variants_by_sample$gene_factor))
  return(variant_frequency)
}


get_tier_frequencies_oncokb <- function(variants_by_sample){
  # makes sure donors can only have a mutation once (no two-knockouts, no two samples)
  by_donor <- variants_by_sample %>% group_by(donor, gene, top_oncokb_level) %>% tally()
  by_donor$n <- 1
  by_donor <- unique(by_donor[,c('donor',   'top_oncokb_level',     'n')])
  
  top_tier_found_for_these_donors <- by_donor$donor[by_donor$top_oncokb_level == "1"]
  
  for(this_tier in c("2", "3B" )){
    
    by_donor$n[by_donor$top_oncokb_level == this_tier & by_donor$donor %in% top_tier_found_for_these_donors] <- 0
    top_tier_found_for_these_donors <- c(top_tier_found_for_these_donors, by_donor$donor[by_donor$top_oncokb_level == this_tier])
  }
  
  variant_frequency <- by_donor %>% group_by(top_oncokb_level) %>% 
    reframe(tally=sum(n), frequency = sum(n) / length(unique(variants_by_sample$donor)))
  
  return(variant_frequency)
}


## replaces get_frequencies from functions to account for "final_status"
get_gene_frequencies_final_status <- function(variants_by_sample){
  # makes sure donors can only have a mutation once (no two-knockouts, no two samples)
  by_donor <- variants_by_sample %>% group_by(donor, gene,  final_status) %>% tally()
  by_donor$n <- 1
  variant_frequency <- by_donor %>% group_by(gene, final_status) %>% 
    reframe(frequency = sum(n) / length(unique(variants_by_sample$donor)))
  
  variant_frequency$gene_factor    <- factor(variant_frequency$gene,  levels = levels(variants_by_sample$gene_factor))
  return(variant_frequency)
}

## replaces get_frequencies from functions to account for "final_status"
get_gene_frequencies_final_status_pathway <- function(variants_by_sample){
  # makes sure donors can only have a mutation once (no two-knockouts, no two samples)
  by_donor <- variants_by_sample %>% group_by(donor, gene, pathway_factor, final_status) %>% tally()
  by_donor$n <- 1
  variant_frequency <- by_donor %>% group_by(gene, pathway_factor, final_status) %>% 
    reframe(frequency = sum(n) / length(unique(variants_by_sample$donor)))
  
  variant_frequency$gene_factor    <- factor(variant_frequency$gene,  levels = levels(variants_by_sample$gene_factor))
  return(variant_frequency)
}

get_gene_frequencies_tissue_type <- function(variants_by_sample ){
  # makes sure donors can only have a mutation once (no two-knockouts, no two samples)
  by_donor <- variants_by_sample %>% group_by(donor, gene, location) %>% tally()
  by_donor$n <- 1
  variant_frequency <- by_donor %>% group_by(gene, location) %>% 
    reframe(frequency = sum(n) / length(unique(variants_by_sample$donor)))
  
  variant_frequency$gene_factor    <- factor(variant_frequency$gene,  levels = levels(variants_by_sample$gene_factor))
  return(variant_frequency)
}

log2_with_minimum <- function(vector, MIN_LOG_ARG=1e-21){
  # """Return log2(x), or log2(min) if x < min; hack to avoid log(0)"""
  
  for(value in vector){
    if (value < MIN_LOG_ARG){
      vector[value] = log2(MIN_LOG_ARG)
    }else{
      vector[value] =  log2(value)
    }
  }
  return(vector)
}

make_empty_action_table <- function(pasted_id){
  var_tmp_df <- cbind.data.frame(
    "donor"= pasted_id,
    "gene"= NA,
    "treatment" = NA,
    "therapy" = NA,
    "rank" = NA
  )
  return(var_tmp_df)
}

make_oncoplot_factor_order <- function(top_variants_by_sample){
  ## takes in a df with columns 'donor' and 'gene'  
  
  gene_tally <- top_variants_by_sample %>% group_by(gene) %>% tally() 
  #refactor 'None' to last
  factor_order <- gene_tally$gene[order(gene_tally$n, decreasing = T)]
  factor_order <- factor_order[factor_order != 'None']
  factor_order <- c(factor_order, 'None')
  gene_tally$gene_factor    <- factor(gene_tally$gene,  levels = factor_order)
  gene_tally <- gene_tally[order(gene_tally$gene_factor),]
  ordered_gene_set <- gene_tally$gene_factor
  number_of_genes = length(ordered_gene_set)
  score_system <- 10^(0:number_of_genes)
  counter_step = number_of_genes
  
  sample_tally <- top_variants_by_sample %>% group_by(donor) %>% tally() 
  sample_tally$order_number <- 0
  
  for (this_gene in ordered_gene_set){
    this_score = score_system[counter_step]
    counter_step = counter_step - 1
    these_samples <- top_variants_by_sample$donor[top_variants_by_sample$gene == this_gene] 
    
    sample_tally$order_number[sample_tally$donor %in% these_samples] <- sample_tally$order_number[sample_tally$donor %in% these_samples] + this_score
    
  }
  
  top_variants_by_sample$donor_factor    <- factor(top_variants_by_sample$donor,  
                                                   levels = sample_tally$donor[order(sample_tally$order_number, decreasing = T)])
  
  top_variants_by_sample$gene_factor    <- factor(top_variants_by_sample$gene,  levels = ordered_gene_set)
  
  return(top_variants_by_sample)
}

make_pvalue_table <- function(immune_melt_IO){
  first=T
  for(this_software in unique(immune_melt_IO$software)){
    tmp_by_software <- immune_melt_IO[immune_melt_IO$software == this_software,]
    for(this_cell_type in unique(tmp_by_software$cell_type)){
      tmp_by_software_by_cell <- tmp_by_software[tmp_by_software$cell_type == this_cell_type,]
      
      test_result <- wilcox.test(tmp_by_software_by_cell$value[tmp_by_software_by_cell$final_status == "No"], 
                                 tmp_by_software_by_cell$value[tmp_by_software_by_cell$final_status == "Yes"], 
                                 alternative = "two.sided")
      pvalue_tmp <- cbind.data.frame(this_software, this_cell_type, signif(test_result$p.value, 2), test_result$statistic)
      if(first){pvalue_table <- pvalue_tmp; first=F}else{
        pvalue_table <- rbind.data.frame(pvalue_table, pvalue_tmp)
      }
    }
  }
  
  names(pvalue_table) <- c("software", "cell_type","pvalue", "W-statistic")
  pvalue_table$sig <- "N.S"
  pvalue_table$sig[pvalue_table$pvalue < 0.05] <- "*"
  pvalue_table$sig[pvalue_table$pvalue < 0.01] <- "**"
  pvalue_table$sig[pvalue_table$pvalue < 0.005] <- "***"
  return(pvalue_table)
}

pull_actions_from_MTB_json <- function(MTB_variants){
  
  first=T
  treatment=therapy=rank=NA
  
  for(req in c(1:length(MTB_variants))){
    
    this_req <- MTB_variants[[req]]
    pasted_id <- concatenate_ids(this_req)
    
    if(length(this_req$mutations) > 0){
      for(mutation in c(1:length(this_req$mutations))){
        
        
        this_mutation <- this_req$mutations[[mutation]]

        if("targetable" %in% names(this_mutation)){
          if(length(this_mutation$targetable) > 0){
            for(action in c(1:length(this_mutation$targetable))){
              actionable_req = T
              
              this_action <- this_mutation$targetable[[action]]
              if("therapy" %in% names(this_action)){ therapy = this_action$therapy} else {therapy = NA}
              if("treatment" %in% names(this_action)){ treatment = this_action$treatment} else {treatment = NA}
              if("rank" %in% names(this_action)){ rank = this_action$rank} else {rank = NA}
              
              if("variant_type" %in% names(this_mutation)){ variant_type = this_mutation$variant_type} else {variant_type = NA}
              if("variant" %in% names(this_mutation)){ variant = this_mutation$variant} else {variant = NA}
              
              var_tmp_df <- cbind.data.frame(
                "BTC_ID"= this_req$ids$BTC_ID,
                "donor"= pasted_id,
                "gene"= this_mutation$gene,
                "variant_type"= variant_type,
                "variant"= variant,
                "treatment" = treatment,
                "therapy" = therapy,
                "rank" = rank
              )
              
              if(first){MTB_df <- var_tmp_df; first=F}else{MTB_df <- rbind.data.frame(MTB_df,var_tmp_df)}
            } 
          } 
        } 
      }
    }
    
  }
  
  return(MTB_df)
}



get_tdp_score <- function(this.sv_file){
  
  #get chromosome sizes
  chrLength_data <- fread('~/Documents/scripts/bitbucket/GRCh38_size.bed')[c(1:22),]
  these.sizes <- chrLength_data$V3[order(chrLength_data$V1)]
  
  dup.chrom.tally <- this.sv_file %>% filter(SV == 'DUP') %>% group_by(CHROM) %>% tally()
  if(nrow(dup.chrom.tally) != 22){
    dup.chrom.tally <-rbind.data.frame(dup.chrom.tally,
                                       cbind.data.frame('CHROM'=chrLength_data$V1[chrLength_data$V1 %ni% dup.chrom.tally$CHROM], 'n'=0)
    )
  }
  dup.chrom.tally <- dup.chrom.tally[order(dup.chrom.tally$CHROM),]
  dup.tally <- this.sv_file %>% filter(SV == 'DUP') 
  
  expected <- nrow(dup.tally) * as.numeric(these.sizes) / sum(these.sizes)
  
  diff <- abs(dup.chrom.tally$n - expected) / nrow(dup.tally)
  
  tdp.score <- (0.71) - sum(diff) ##published formula adjustment
  
  return(tdp.score)
}

pull_expression <- function(sample_list, gene_list, basedir = "/Volumes/pcsi/pipeline/hg38_production/" ){
  names(sample_list) <- paste0('V',c(1:length(names(sample_list))))
  
  first=T
  for(sample_id in sample_list$V2){
    
    donor_id = unique(sample_list$V1[sample_list$V2 == sample_id])
    
    file_path=unique(paste0(basedir, donor_id,"/", sample_id,"/rna/star/2.7.4a/stringtie/2.0.6/",sample_id, "_stringtie_abundance.txt"))
    
    if(file.exists(file_path)){
      cat("found rna for ",sample_id, "\n")
      gene_file <- fread(file_path)[,c("Gene Name", "TPM")]
      names(gene_file)[1] <- c('gene_name')
      
      gene_group <- gene_file %>%
        group_by(gene_name) %>%
        summarise(TPM = sum(TPM))
      
      gene_filtered <- gene_group %>% filter( gene_name %in% gene_list)
      gene_filtered$tpm_adj <- gene_filtered$TPM / sum(gene_filtered$TPM) * 1000000
      
      names(gene_filtered)[3] <- c(sample_id)
      
      if(first){gene_df <- as.data.frame(gene_list); joined_file <- left_join(gene_df, gene_filtered[,c(1,3)], by=c("gene_list"="gene_name")); first=F}else{
        
        joined_file <- left_join(joined_file, gene_filtered[,c(1,3)], by=c("gene_list"="gene_name"))
        
      }
      
    } else {
      cat("rna not found for ",sample_id, "\n")
      
    }
    
  }

  return(joined_file)
}

pull_mutations_from_MTB_json <- function(MTB_variants){
  first=T
  for(req in c(1:length(MTB_variants))){
    
    this_req <- MTB_variants[[req]]
    cat(req," ")
    
    pasted_id <- concatenate_ids(this_req)
    
    if(length(this_req$mutations) > 0){
      for(mutation in c(1:length(this_req$mutations))){
        this_mutation <- this_req$mutations[[mutation]]
        
        var_tmp_df <- cbind.data.frame(
          "BTC_ID"= this_req$ids$BTC_ID,
          "donor"= pasted_id,
          "gene"= this_mutation$gene
        )
        cat(mutation,"\n")
        print(var_tmp_df)
        if(first){MTB_df <- var_tmp_df; first=F}else{MTB_df <- rbind.data.frame(MTB_df,var_tmp_df)}
        
      }
    }else{
      var_tmp_df <- cbind.data.frame(
        "BTC_ID"= this_req$ids$BTC_ID,
        "donor"= pasted_id,
        "gene"= NA
      )
      cat(mutation,"\n")
      print(var_tmp_df)
      if(first){MTB_df <- var_tmp_df; first=F}else{MTB_df <- rbind.data.frame(MTB_df,var_tmp_df)}
      
    }
    
  }
  
  return(MTB_df)
}

pull_rip_from_MTB_json <- function(MTB_variants){
  first=T
  for(req in c(1:length(MTB_variants))){
    
    this_req <- MTB_variants[[req]]
    
    pasted_id <- concatenate_ids(this_req)
    rip = "Alive at MTB"
    if("hx" %in% names(this_req)){
      for(history in c(1:length(this_req$hx))){
        this_history = this_req$hx[[history]]
        if("RIP" %in% names(this_history)){
          rip = "Succumbed by MTB" 
        }
      }
    }
    var_tmp_df <- cbind.data.frame(
      "BTC_ID"= this_req$ids$BTC_ID,
      "pasted_id"= pasted_id,
      "rip"= rip
    )
    if(first){rip_df <- var_tmp_df; first=F}else{rip_df <- rbind.data.frame(rip_df,var_tmp_df)}
    
  }
  return(rip_df)
}

pull_variants_from_CAP_json <- function(CAP_variants){
  first=T
  for(req in c(1:length(CAP_variants))){
    this_req <- CAP_variants[[req]]
    
    if(this_req$QC == "PASS" & length(this_req$mutations) > 0){
      for(mutation in c(1:length(this_req$mutations))){
        this_mutation <- this_req$mutations[[mutation]]
        
        var_tmp_df <- cbind.data.frame(
          "donor"= this_req$study_id,
          "gene"= this_mutation$gene,
          "variant"= this_mutation$variant,
          "oncokb_level"= this_mutation$oncokb_level
        )
        
        if(first){CAP_df <- var_tmp_df;first=F}else{CAP_df <- rbind.data.frame(CAP_df,var_tmp_df)}
        
      }
    }
  }
  return(CAP_df)
}

check_specific_substitutions <- function(these_variants, this_alteration, actionable_variants_list, this_level,this_indication, this_drug, this_gene){
  for(variant_row in c(1:nrow(these_variants))){
    
    if(!is.na(these_variants$aa_context[variant_row])){
      
      # split possible transcripts by "|"
      these_transcripts <- strsplit(these_variants$aa_context[variant_row],"|",fixed=TRUE)
      for(this_transcript in these_transcripts[[1]]){
        # split transcript ID from substitution
        this_substitution <- strsplit(this_transcript,":p.",fixed=TRUE)[[1]][2]
        
        if(grepl(this_alteration, this_substitution, fixed = TRUE) ){
          actionable_variants_list <- list.append(actionable_variants_list, 
                                                  list(these_variants$donor[variant_row], this_alteration, this_level,this_indication, this_drug, this_gene, this_substitution))
        }
      }
    }
  }
  return(actionable_variants_list)
}


get_first_aa <- function(these_variants){
  these_aas <- c()
  for(variant_row in c(1:nrow(these_variants))){
    
    if(!is.na(these_variants$aa_context[variant_row])){
      
      # split possible transcripts by "|"
      these_transcripts <- strsplit(these_variants$aa_context[variant_row],"|",fixed=TRUE)

      # split transcript ID from substitution
      this_substitution <- strsplit(these_transcripts[[1]],":p.",fixed=TRUE)[[1]][2]
      these_aas <- c(these_aas, this_substitution)
      
    } 
  }
  return(these_aas)
}


get_oncokb_actionable_variants <- function( somatic_variants_oncoKB_filtered, oncokb_genes, oncogenic_variants, TRUNCATING_MUTATIONS, INCLUDES_BILIARY, cohort=NULL, summary=NULL, actionable_fusions=NULL, other_biomarkers=NULL){
  require(rlist)
  
  actionable_variants_list <- list()
  for(row in c(1:nrow(oncokb_genes))){
    
    this_row <- oncokb_genes[row,]
    this_gene = this_row$Gene
    cat("\nchecking row ", row,": ",this_gene, " ",this_row$Alterations, "\n")
    these_alterations <- strsplit(this_row$Alterations,", ")
    these_variants <- somatic_variants_oncoKB_filtered[somatic_variants_oncoKB_filtered$gene == this_gene,]
    this_level <- this_row$Level
    this_indication <- this_row$`Cancer Types`
    if(this_indication %ni% INCLUDES_BILIARY){
      this_level <- '3.5'
    }
    this_drug <- this_row$`Drugs (for therapeutic implications only)`
    
    if (this_row$Alterations  == "Microsatellite Instability-High"){
      #not currently in variant list, will need another approach
      
      if(!is.null(cohort) & !is.null(summary)){
        for(this_donor in cohort$BTCID){
          if(this_donor %in% summary$donor  ){
            if(unique(summary$mmr_bin[summary$donor == this_donor]) == "MSI"){
            actionable_variants_list <- list.append(actionable_variants_list, 
                                                    list(this_donor, "MSI-H", this_level,this_indication, this_drug, "MSI", "MSI"))
            
          }
          } else {
            cat(this_donor,' not in summary file\n')
          }
        }
      } else {
        cat("no summary file, skipping MSI\n")
      }
    } else if(this_row$Alterations == "Tumor Mutational Burden-High"){
      if(!is.null(cohort) & !is.null(summary)){
        
        for(this_donor in cohort$BTCID){
          if(this_donor %in% summary$donor  ){
             if( (mean(summary$snv_count[summary$donor == this_donor], na.rm = T) + mean(summary$indel_count[summary$donor == this_donor], na.rm = T)) / 3000 > 10){
               actionable_variants_list <- list.append(actionable_variants_list, 
                                                       list(this_donor, "TMB-H", this_level,this_indication, this_drug, "TMB", "TMB"))
               
             }
          } else {
            cat(this_donor,' not in summary file\n')
          }
        }
      }else {
        cat("no summary file, skipping TMB\n")
      }
      
    } else if(nrow(these_variants) == 0){ cat("no variants found in this gene\n") 
    } else if(length(these_alterations[[1]]) > 1){
      #cat("checking multiple possible alterations:")
      for(this_alteration in these_alterations[[1]]){
        #cat(" checking for ", this_alteration, "\n")
        actionable_variants_list <- check_specific_substitutions(these_variants, this_alteration, actionable_variants_list, this_level,this_indication, this_drug, this_gene)
      }
    } else {
      
      this_alteration <-  these_alterations[[1]]
      
      if(this_alteration == "Amplification"){
        
        for(variant_row in c(1:nrow(these_variants))){
          if(these_variants$mutation_type[variant_row] == "strong amplification"){
            cat(these_variants$donor[variant_row], "\n")
            actionable_variants_list <- list.append(actionable_variants_list, 
                                                    list(these_variants$donor[variant_row], this_alteration, this_level,this_indication, this_drug, this_gene, these_variants$mutation_type[variant_row]))
            
          }
        }
        
        if(!is.null(other_biomarkers)){
          ## CHECK AMPs from other sources
          these_other_biomarkers <- other_biomarkers[other_biomarkers$gene == this_gene,]
          
          if(nrow(these_other_biomarkers) > 0){
            for(variant_row in c(1:nrow(these_other_biomarkers))){
              if(these_other_biomarkers$alteration_class[variant_row] == "Amplification"){
                cat(these_other_biomarkers$btc_id[variant_row], "\n")
                actionable_variants_list <- list.append(actionable_variants_list, 
                                                        list(these_other_biomarkers$btc_id[variant_row], this_alteration, this_level, this_indication, this_drug, this_gene, these_other_biomarkers$source[variant_row]))
                
              }
            }
          }
        } else {
          cat("no other biomarkers given, skipping...\n")
        }
        
      } else if (this_alteration == "Fusions"){
        ## CHECK RNA for fusions
        if(!is.null(actionable_fusions)){
          these_actionable_fusions <- actionable_fusions[actionable_fusions$gene == this_gene,]
        
          if(nrow(these_actionable_fusions) > 0){
            for(variant_row in c(1:nrow(these_actionable_fusions))){
                ## MAKE SURE THIS FUSION MATCHES RNA FUSION
                
                #cat(actionable_fusions$donor[variant_row],actionable_fusions$gene[variant_row], "\n")
                actionable_variants_list <- list.append(actionable_variants_list, 
                                                        list(these_actionable_fusions$donor[variant_row], this_alteration, this_level, this_indication, this_drug, this_gene, these_actionable_fusions$gene[variant_row]))
                
              
            }
          }
        } else{
          cat("no fusion file given, skipping fusions\n")
        }
        
      } else if (this_alteration == "Truncating Mutations"){
        for(variant_row in c(1:nrow(these_variants))){
          cat(these_variants$donor[variant_row]," ",these_variants$gene[variant_row]," ",these_variants$mutation_type[variant_row], "\n")
          if(these_variants$mutation_type[variant_row] %in% TRUNCATING_MUTATIONS){
            actionable_variants_list <- list.append(actionable_variants_list, 
                                                    list(these_variants$donor[variant_row], this_alteration, this_level,this_indication, this_drug, this_gene, these_variants$mutation_type[variant_row]))
          }
        }
        
      } else  if (this_alteration == "Oncogenic Mutations"){
        these_oncogenic_variants <- oncogenic_variants$variant[oncogenic_variants$gene == this_gene]
        for(variant_row in c(1:nrow(these_variants))){
          cat(these_variants$donor[variant_row]," ",these_variants$gene[variant_row]," ",these_variants$mutation_type[variant_row], "\n")
          
          for(this_oncogenic_variant in these_oncogenic_variants){
            if(this_oncogenic_variant %in% c("Deletion", "Amplification")){
              if(these_variants$mutation_type[variant_row] %in% c("strong amplification", "homozygous deletion")){
                #cat(these_variants$donor[variant_row], "\n")
                actionable_variants_list <- list.append(actionable_variants_list, 
                                                        list(these_variants$donor[variant_row], this_alteration, this_level,this_indication, this_drug, this_gene, these_variants$mutation_type[variant_row]))
                
              }
            }else {
              actionable_variants_list <- check_specific_substitutions(these_variants[variant_row,], this_oncogenic_variant, actionable_variants_list, this_level,this_indication, this_drug, this_gene)
            }
          }
        }
        
      } else {
        # else assumes specific substitutions
        actionable_variants_list <- check_specific_substitutions(these_variants, this_alteration, actionable_variants_list, this_level,this_indication, this_drug, this_gene)
      }
    }
  }
  
  actionable_variants_df <- do.call(rbind.data.frame, actionable_variants_list)
  names(actionable_variants_df) <- c("donor", "alteration_class", "oncokb_level","indication", "treatment", "gene", "alteration_type")
  
  return(actionable_variants_df)
}



get_actionable_mutations <- function(psc_actionable, actionable_list, mutation_type.actionable_type ){
  first=T
  for(row in c(1:nrow(actionable_list))){
    this_row = actionable_list[row,]
    this_gene = this_row$gene
    these_mutation_types = 
      mutation_type.actionable_type$mutation_type[mutation_type.actionable_type$actionable_type == this_row$actionable_type]
    
    if(first){
      actionable_variants_df <- 
        psc_actionable[psc_actionable$gene == this_gene & 
                         psc_actionable$mutation_type %in% these_mutation_types,]
      first=F
    }else{
      actionable_variants_df <- 
        rbind.data.frame(psc_actionable[psc_actionable$gene == this_gene & 
                                          psc_actionable$mutation_type %in% these_mutation_types,],
                         actionable_variants_df)
      
    }
    
  }
  actionable_variants <- left_join(actionable_variants_df, mutation_type.actionable_type)
  
  return(actionable_variants)
}

add_sv_clusters_to_summary <- function(chrLength_data, summary_path="/Volumes/pcsi/users/fbeaudry/BTC.summary.csv"){
  summary_raw <- fread(summary_path)
  summary_raw$assay <- '80X'
  summary_raw$assay[summary_raw$tumour_coverage < 60] <- '40X'
  
  for(this_sample_id in summary_raw$tumour){ 
    this_donor_id <- summary_raw$donor[summary_raw$tumour == this_sample_id]
    file_path <- paste0('/Volumes/pcsi/pipeline/hg38/',this_donor_id,'/',this_sample_id,'/wgs/bwa/0.7.17/mergeSV/',this_sample_id,'_finalannotatedSV.txt')
    
    if(file.exists(file_path)){
      cat(this_sample_id, "\n")
      finalannotatedSV <- fread(file_path)
      
      chrLength <- chrLength_data
      chrLength$tally <- 0
      
      # tally variants per chromosome
      for (this_sv_index in c(1:nrow(finalannotatedSV))) {
        
        chrLength$tally[chrLength$V1 == finalannotatedSV$chrom1[this_sv_index]] <- 
          chrLength$tally[chrLength$V1 == finalannotatedSV$chrom1[this_sv_index]] + 1
        
      }
      
      #normalize by mb
      chrLength$sv_per_mb = chrLength$tally / (chrLength$V3 / 1000000)
      
      #sort by number of sv / mb
      rankedSvPerChr = chrLength[order(chrLength$sv_per_mb),]
      
      ## from Waddell et al 2015:
      ## "chromosomes with a high structural variant mutation rate per Mb 
      ## exceeded 5 times the length of the interquartile range 
      ## from the 75th percentile of the chromosome counts for each patient"
      
      svPerMinimum = 5 * ((rankedSvPerChr$sv_per_mb[17] + rankedSvPerChr$sv_per_mb[18]) / 2)
      
      ## adjusted quartiles to instead represent chromoplexy 
      extreme_svPerMinimum = 5 * ((rankedSvPerChr$sv_per_mb[20] + rankedSvPerChr$sv_per_mb[21]) / 2)
      
      summary_raw$wadell_factor[summary_raw$tumour == this_sample_id] <- max(rankedSvPerChr$sv_per_mb) - svPerMinimum
      summary_raw$sv_cluster_factor[summary_raw$tumour == this_sample_id] <- max(rankedSvPerChr$sv_per_mb) - extreme_svPerMinimum
    }
  }
  
  return(summary_raw)
}


predict_TSP_with_confidence <- function(tsp_classifier, tpm_matrix, confidence_cutoff = 0.75) {
  
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

