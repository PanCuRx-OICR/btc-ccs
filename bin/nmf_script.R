# nmf_script.R
# mod. from Duhan's script, forked on May 22, 2025

library(NMF)
library(data.table)
library(dplyr)

core.dir = '/.mounts/labs/PCSI/'
#core.dir = '/Volumes/pcsi/'

outdir = file.path(core.dir,'users/fbeaudry/nmf_analysis/')
MIN_READ = 0

#### filter for analysis ####

gene_set <- fread(file.path(outdir,'/fgs.txt'))
gene_list <- gene_set$ensid 

joined_file_filt <- fread(file.path(outdir,"/rna.read.ff.txt"))

protein_genes <- gene_set$gene_name[gene_set$biotype=="protein_coding" & gene_set$Length >= 500 & gene_set$Chr %in% c(1:23)]
joined_file_protein <- joined_file_filt %>% filter(gene_name %in% protein_genes)

genes <- joined_file_protein$gene_name
btc_tpm <- as.matrix(joined_file_protein[,-1])
rownames(btc_tpm) <- genes

min_samples <- round(0.50*ncol(btc_tpm))

btc_tpm_filtered <- btc_tpm[rowSums(btc_tpm > MIN_READ) >= min_samples,]

btc_tpm_filtered <- log(btc_tpm_filtered + 1)

# Preprocess
for(N_GENES in c(2000, 4000, 6000, 8000, 10000, 12000, 14000)){ 
  
  cat(N_GENES,'\t')
  
  # get the variances
  gene_variances <- apply(btc_tpm_filtered, 1, var)
  
  # select the top genes for downstream analysis
  top_genes <- order(gene_variances, decreasing = TRUE)[1:N_GENES]
  
  # save the highly variable genes for NMF
  filtered_genes_cts <- na.omit(btc_tpm_filtered[top_genes, ])
  filtered_genes_cts <- as.matrix(filtered_genes_cts)
  saveRDS(filtered_genes_cts, file = paste0(outdir,"nmf.",N_GENES,"g.rds"))
  
  #### run NMF ####
  
  for(component_n in c(2:20)){
    cat(component_n,'\t')
  
    nmf_res <- nmf(x = filtered_genes_cts,
    	       rank = component_n,
    	       method = 'lee',
    	       seed = 1234,
    	       nrun = 200, 
    	       .options = "t")
    
    H_matrix = coef(nmf_res)
    W_matrix = basis(nmf_res)
    
    # save all the outputs
    
    saveRDS(nmf_res, file = paste0(outdir,"nmf.k",component_n,".",N_GENES,"g.rds"))
    write.csv(H_matrix, paste0(outdir,"h_mat.k",component_n,".",N_GENES,"g.csv"))
    write.csv(W_matrix, paste0(outdir,"w_mat.k",component_n,".",N_GENES,"g.csv"))
  }

}





