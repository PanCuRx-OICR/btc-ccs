
#BiocManager::install("CNTools")

suppressPackageStartupMessages({
  
require(mclust)
require(data.table)
require(ggplot2)
require(dplyr)
require(tidyr)
require(CNTools)
  
})

add_abs_cn <- function(this.gene, this.cnv.tpm,  gene.cns){
  this.cn <- data.frame('acn'=as.numeric(t(gene.cns[gene.cns$genename == this.gene, ])[-c(1:5),]))
  this.cn$tumour <- names(gene.cns)[-c(1:5)]
  
  suppressMessages(
    this.cn.tpm <- inner_join( this.cnv.tpm, this.cn,by=c('tumour'='tumour'))
  )
  return(this.cn.tpm)
}

add.bin.n <- function(this.cn.tpm, tpm.glm.coef, call_type){
  this.call.n <- this.cn.tpm %>% group_by(call) %>% tally()
  if(!is.na(tpm.glm.coef$bin.thresh)){
    tpm.glm.coef$bin.n <- this.call.n$n[this.call.n$call == call_type]
  } else {
    tpm.glm.coef$bin.n <- NA
  }
  return(tpm.glm.coef)
}

add_bin_threshold <- function(this.cn.tpm, tpm.glm.coef, call_type){
  this.bin <- this.cn.tpm %>% group_by( cn) %>% summarise(mean_bin = mean(call)) %>% filter(mean_bin > 0.5)
  if(nrow(this.bin) > 0){
    if(call_type == 'DEL'){
      
      this.bin <- this.bin %>% filter( cn == max(cn))
    } else if(call_type == 'AMP'){
      
      this.bin <- this.bin %>% filter( cn == min(cn))
    }
    tpm.glm.coef$bin.thresh <- this.bin$cn
  } else {
    tpm.glm.coef$bin.thresh <- NA
  }
  return(tpm.glm.coef)
}

add_tpm_clusters <- function(this.cn.tpm, LOG_TPM){
  
  suppressMessages(
    if(LOG_TPM){
      fit_mclust <- Mclust(log(this.cn.tpm$tpm + 1),G=2,verbose = F)
    } else {
      fit_mclust <- Mclust(this.cn.tpm$tpm,G=2,verbose = F)
    }
  )
  
  this.cn.tpm$class <- fit_mclust$classification
  
  this.cn.tpm$cn <- round(this.cn.tpm$acn)
  return(this.cn.tpm)
}

adjust_for_ploidy <- function(this.cn.tpm, ploidies){
  this.cn.tpm.ploidies <- left_join(this.cn.tpm, ploidies, by=c('tumour'='tumour'))
  this.cn.tpm.ploidies$cn[this.cn.tpm.ploidies$ploidy.bin == 'WGD'] <- round(this.cn.tpm.ploidies$acn[this.cn.tpm.ploidies$ploidy.bin == 'WGD'] / 2)
  this.cn.tpm.ploidies <- this.cn.tpm.ploidies %>% dplyr::select(-ploidy.bin)
  return(this.cn.tpm.ploidies)
}

get.sv.window <- function(this.gene, gene_set, WINDOW_SIZE, unclustered.sv){
  this.gene.chr <- gene_set$Chr[gene_set$gene_name == this.gene]
  window.start <- gene_set$Start[gene_set$gene_name == this.gene] - WINDOW_SIZE
  window.end <- gene_set$End[gene_set$gene_name == this.gene] + WINDOW_SIZE
  
  rearranged.tumours <- unclustered.sv[unclustered.sv$CHROM == this.gene.chr & window.start < unclustered.sv$POS & unclustered.sv$POS < window.end, ]
  
  rearranged.tumours$GENE[is.na(rearranged.tumours$GENE)] <- this.gene
  rearranged.tumours$SV[is.na(rearranged.tumours$SV)]  <- 'BND'
  
  return(rearranged.tumours)
}

make_cn_tpm_plots <- function(this.gene, this.cn.tpm, root_path){
  
  png(paste0( root_path,"/",this.gene,".bin.png") ,width = 800, height = 800) #, bg = "transparent"
  
  print(
    
    ggplot(this.cn.tpm, aes(x=factor(cn), y=log(tpm)) )+ 
      geom_boxplot(outlier.shape = NA) + 
      geom_jitter(aes(color=factor(call), shape=factor(class)), size=2, width = 0.2) + 
      guides()+
      scale_color_manual(values=c('cornflowerblue','forestgreen')) +
      scale_shape_manual(values=c(1,16))+
      theme_bw(base_size=22) + 
      labs(x='CN', y='TPM', title=this.gene, color='Call', shape='TPM.bin') +
      theme(panel.grid = element_blank()) 
    
  )
  
  dev.off()
  
  svg(paste0( root_path,"/",this.gene,".bin.svg") ,width = 4, height = 6) #, bg = "transparent"
  
  print(
    
    ggplot(this.cn.tpm, aes(x=factor(call), y=log(tpm)) )+ 
      geom_boxplot(outlier.shape = NA) + 
      geom_jitter(aes(shape=factor(class)), size=2, width = 0.2) + 
      guides()+
      #scale_color_manual(values=c('cornflowerblue','forestgreen')) +
      scale_shape_manual(values=c(1,16))+
      theme_bw(base_size=22) + 
      labs(x='CN', y='TPM', title=this.gene, shape='Class') +
      theme(panel.grid = element_blank()) 
    
  )
  
  dev.off()
  
  
  png(paste0( root_path,"/",this.gene,".png") ,width = 600, height = 600) #, bg = "transparent"
  
  print(
    ggplot(this.cn.tpm, aes(x=cnv, y=log(tpm))) + 
      
      geom_point() + 
      geom_smooth(method = lm, formula = y ~ x)+
      theme_bw(base_size=22) + 
      labs(x='CNV log2(tumor/normal)', y='TPM', title=this.gene) +
      theme(panel.grid = element_blank()) 
  )
  
  dev.off()
}

make_cnv_tpm_table <- function(this.gene, rna, gene.cn){
  
  this.expression <- data.frame('tpm'=as.numeric(t(rna[rna$gene_list == this.gene, ])[-c(1,2),]))
  this.expression$tumour <- names(rna)[-c(1:2)]
  
  this.cnv <- data.frame('cnv'=as.numeric(t(gene.cn[gene.cn$genename == this.gene, ])[-c(1:5),]))
  this.cnv$tumour <- names(gene.cn)[-c(1:5)]
  
  this.cnv.tpm <- inner_join(this.cnv, this.expression, by=c('tumour'='tumour'))
  return(this.cnv.tpm)
}

make_empty_tpm_table <- function(tpm.glm.coef){
  tpm.glm.coef$cluster.p <- NA
  tpm.glm.coef$bin.thresh <- NA
  tpm.glm.coef$bin.n <- NA
  tpm.glm.coef$bin.p <- NA
  tpm.glm.coef$bin.w <- NA
  tpm.glm.coef$bin.est <- NA
  tpm.glm.coef$bin.test.p <- NA
  
  return(tpm.glm.coef)
}

match_acn_to_clusters <- function(this.cn.tpm, call ){
  class.means <- this.cn.tpm %>% group_by(class) %>% summarise(class.mean = mean(cnv))
  
  if(call == 'DEL'){
    this_class <- class.means %>% filter( class.mean == min(class.mean))
  } else if(call == 'AMP'){
    this_class <- class.means %>% filter( class.mean == max(class.mean))
  }
  
  this.cn.tpm$call <- 0
  this.cn.tpm$call[this.cn.tpm$class == this_class$class] <- 1
  
  
  return(this.cn.tpm)
}

run_gistic_filter <- function(these.genes, rna, gene.cn, 
                              gene.cns, basedir, call_type, 
                              LOG_TPM,  linear.p = NULL, ploidies = NULL, 
                              linear.test.type = 'sc'){
  
  if(is.null(linear.p)){
    #bonferroni
    linear.p = (0.05)/length(these.genes)
  }
  
  first.loop=T
  first.gene=T
  
  for(this.gene in these.genes){
    
    this.cnv.tpm <- make_cnv_tpm_table(this.gene, rna, gene.cn)
    
    active.expression <- this.cnv.tpm %>% filter(tpm > 0)
    cn.called.loci <- this.cnv.tpm %>% filter(cnv > 0)
    
    if(nrow(active.expression) > 1 & nrow(cn.called.loci) > 1 ){ 
      
      if(linear.test.type == 'lm'){
        tpm.glm.coef <- run_linear_model(this.cnv.tpm, LOG_TPM)
      } else if(linear.test.type == 'sc'){
        tpm.glm.coef <- run_rank_correlation(this.cnv.tpm)
      }
      
      tpm.glm.coef$gene <- this.gene
      tpm.glm.coef$call <- call_type
      
      if( tpm.glm.coef$`Pr(>|t|)` < linear.p ){
        
        this.cn.tpm <- add_abs_cn(this.gene, this.cnv.tpm,  gene.cns)
        this.cn.tpm <- add_tpm_clusters(this.cn.tpm, LOG_TPM)
        
        if(!is.null(ploidies)){
          
          this.cn.tpm <- adjust_for_ploidy(this.cn.tpm, ploidies)
          
        }
        
        this.cn.tpm <- match_acn_to_clusters(this.cn.tpm, call = call_type)
        
        tpm.glm.coef <- run_wilcox_cluster_test(this.cn.tpm, tpm.glm.coef)
        tpm.glm.coef <- add_bin_threshold(this.cn.tpm, tpm.glm.coef, call_type)
        
        this.cn.tpm$call <- 'WT'
        
        if(call_type == 'DEL'){
          
          #deletions are less than or equal to threshold
          this.cn.tpm$call[this.cn.tpm$cn <= tpm.glm.coef$bin.thresh] <- call_type
          
        } else if (call_type == 'AMP') {
          
          this.cn.tpm$call[this.cn.tpm$cn >= tpm.glm.coef$bin.thresh] <- call_type
          
        }
        
        this.cn.tpm$call <- factor(this.cn.tpm$call, levels=c('DEL','WT','AMP' ))
        
        tpm.glm.coef <- add.bin.n(this.cn.tpm, tpm.glm.coef, call_type)
        
        
        if(all(!is.na(tpm.glm.coef$bin.n) & tpm.glm.coef$bin.n > 1 & tpm.glm.coef$bin.n < (nrow(this.cn.tpm)-1)  )){
          
          
          if(LOG_TPM){
            w.bin.tst.res <- wilcox.test(log(this.cn.tpm$tpm[this.cn.tpm$call == 'WT'] + 1), log(this.cn.tpm$tpm[this.cn.tpm$call != 'WT'] + 1), conf.int = T)
          }else{
            w.bin.tst.res <- wilcox.test(this.cn.tpm$tpm[this.cn.tpm$call == 'WT'], this.cn.tpm$tpm[this.cn.tpm$call != 'WT'], conf.int = T)
          }
          
          tpm.glm.coef$bin.p <- w.bin.tst.res$p.value
          tpm.glm.coef$bin.w <- w.bin.tst.res$statistic
          tpm.glm.coef$bin.est <- w.bin.tst.res$estimate
          
          #confirm that binning went OK
          w.bin.chk.res <- wilcox.test(this.cn.tpm$cn[this.cn.tpm$call == 'WT'], this.cn.tpm$cn[this.cn.tpm$call != 'WT'])
          
          tpm.glm.coef$bin.test.p <- w.bin.chk.res$p.value
          
          this.call <- this.cn.tpm %>% dplyr::select(tumour,call)
          names(this.call)[2] <- this.gene
          if(first.gene){all.calls <- this.call; first.gene=F}else{all.calls <- full_join(all.calls,this.call,by=c('tumour'='tumour'))}
          
          plot_path = file.path(basedir,"/lm/",call_type)
          
          if (!file.exists(plot_path)) {
            # Create the directory if it doesn't exist
            dir.create(plot_path, recursive = TRUE)
            cat("Directory created:", plot_path, "\n")
          } 
          
          cat('plotting ', this.gene,'\n')
          make_cn_tpm_plots(this.gene, this.cn.tpm, root_path=plot_path)
          
          
        } else {
          
          tpm.glm.coef$bin.p <- NA
          tpm.glm.coef$bin.w <- NA
          tpm.glm.coef$bin.est <- NA
          tpm.glm.coef$bin.test.p <- NA
          
        }
        
        
      }else{
        #if CNV not significantly associated with TPM
        tpm.glm.coef <- make_empty_tpm_table(tpm.glm.coef)
      }
      
      if(first.loop){del.tpm <- tpm.glm.coef; first.loop=F} else {del.tpm <- rbind.data.frame(del.tpm, tpm.glm.coef)}
    }
    
  }
  
  write.table(
    del.tpm,
    file = paste0(basedir, "/tpm.cn.",call_type,".txt"),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
  write.table(
    unique(all.calls),
    file = paste0(basedir, "/",call_type,".calls.txt"),
    row.names = FALSE, quote = FALSE, sep = "\t"
  )
}

run_linear_model <- function(this.cnv.tpm, LOG_TPM){
  if(LOG_TPM){
    tpm.glm <- lm(log(as.numeric(tpm) + 1) ~ as.numeric(cnv), data=this.cnv.tpm)
  } else {
    tpm.glm <- lm(as.numeric(tpm) ~ as.numeric(cnv), data=this.cnv.tpm)
  }
  
  tpm.glm.summary <- summary(tpm.glm)
  if(nrow(tpm.glm.summary$coefficients) > 1){
    tpm.glm.coef <- as.data.frame(t(tpm.glm.summary$coefficients[2,]))[-2]
  } else {
    tpm.glm.coef <-cbind.data.frame(
      'Estimate'   =NA,
   #   'Std. Error' =NA,
      't value'    =NA,
      'Pr(>|t|)'   =NA
    )
  }
  return(tpm.glm.coef)
}

run_rank_correlation <- function(this.cnv.tpm){
  #Spearman's rank correlation
  
  rank_cnv <- rank(this.cnv.tpm$cnv)
  rank_tpm <- rank(this.cnv.tpm$tpm)
  
  spearman_corr <- cor.test(rank_cnv, rank_tpm, method = "pearson")
  
  spearman.coef <-cbind.data.frame(
    'Estimate'   =spearman_corr$estimate,
    't value'    =spearman_corr$statistic,
    'Pr(>|t|)'   =spearman_corr$p.value
  )
  return(spearman.coef)
}

run_wilcox_cluster_test <- function(this.cn.tpm, tpm.glm.coef){
  
  w.tst.res <- wilcox.test(this.cn.tpm$cn[this.cn.tpm$class == 1], this.cn.tpm$cn[this.cn.tpm$class == 2])
  tpm.glm.coef$cluster.p <- w.tst.res$p.value
  #tpm.glm.coef$cluster.w <- w.tst.res$statistic
  return(tpm.glm.coef)
}

split_and_sum <- function(this.value) {
  parts <- strsplit(as.character(this.value), "\\.")[[1]]
  sum(as.numeric(parts))
}









