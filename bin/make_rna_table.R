```{r}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

basedir = '~/Documents/scripts/github/btc-ccs/'
secdir = '~/Documents/scripts/github/evolution/'
source( file.path(basedir,"/bin/BTC.functions.R"))

#### from htseq ####

hpc.dir = '/Volumes/pcsi/'

outdir = file.path(hpc.dir,'users/fbeaudry/nmf_analysis/')
pipepeline_base_dir = file.path(hpc.dir,"/pipeline/hg38_production/")

gene_set <- fread(file.path(basedir,'/source_data/fgs.txt'))
gene_set$kb <- gene_set$Length / 1000

## no filtering of genes to ensure better portability to signature classifiers
full_gene_list <- gene_set$ensid 

lm_set <- fread(file.path(basedir,'/source_data/LM22.txt'))

narrow_gene.df <- gene_set %>% filter( ( biotype=="protein_coding" & Length >= 500 & Chr %in% c(1:23) ) | 
                                  gene_name %in% lm_set$`Gene symbol`) %>% dplyr::select(ensid, gene_name, kb) 

sample_list <- fread(file.path(hpc.dir,'/users/fbeaudry/every.sample.list.txt'), header=F)
sample_list$V1 <- gsub("_","",sample_list$V1)
sample_list[sample_list == ''] <- NA
sample_list <- sample_list %>% filter(!is.na(V1) & !is.na(V3))

pull.rna = F
if(pull.rna){
  first=T
  for(sample_id in sample_list$V3){
    
    donor_id= sample_list$V1[sample_list$V3 == sample_id]
    
    file_path=paste0(pipepeline_base_dir, donor_id,"/", sample_id,"/rna/star/2.7.4a/htseq/hg38_ens100/",sample_id, "_htseq_count_all.txt")
    
    if(file.exists(file_path)){
      cat("found rna for ",sample_id, "\n")
      gene_file <- fread(file_path) 
      
      names(gene_file) <- c('gene_name', sample_id)
      
      gene_filtered <- gene_file %>% filter( gene_name %in% full_gene_list)
      
      if(first){
        first=F
        
        #gene_df <- data.frame('gene_list'=full_gene_list)
        these.counts <- left_join(gene_set, gene_filtered, by=c("ensid"="gene_name"))
        names(these.counts)[10] <- 'sample'
        joined_file <- these.counts %>% group_by(gene_name) %>% 
          summarise(counts = sum(sample)) %>%
          dplyr::select(all_of(c('gene_name','counts')))
        names(joined_file)[2] <- sample_id
        
        tpm_data <- left_join(narrow_gene.df, gene_filtered, by=c("ensid"="gene_name"))
        tpm_data[,4] <- tpm_data[,4]/tpm_data$kb
        scaling_factors <- colSums(tpm_data[,4])
        tpm_data[,4] <- (tpm_data[,4]/ scaling_factors) *  1e6
        
        names(tpm_data)[4] <- 'sample'
        tpm_data <- tpm_data %>% group_by(gene_name) %>% summarise(tpm = sum(sample))
        names(tpm_data)[2] <- sample_id
        
        all_rna_raw <- tpm_data
        
        
      }else{
        
        these.counts <- left_join(gene_set, gene_filtered, by=c("ensid"="gene_name"))
        names(these.counts)[10] <- 'sample'
        these.counts <- these.counts %>% group_by(gene_name) %>% 
          summarise(counts = sum(sample)) %>%
          dplyr::select(all_of(c('gene_name','counts')))
        names(these.counts)[2] <- sample_id
        
        joined_file <- left_join(joined_file, these.counts, by=c("gene_name"="gene_name"))
        
        tpm_data <- left_join(narrow_gene.df, gene_filtered, by=c("ensid"="gene_name"))
        tpm_data[,4] <- tpm_data[,4]/tpm_data$kb
        scaling_factors <- colSums(tpm_data[,4])
        tpm_data[,4] <- (tpm_data[,4]/ scaling_factors) *  1e6
        
        names(tpm_data)[4] <- 'sample'
        tpm_data <- tpm_data %>% group_by(gene_name) %>% summarise(tpm = sum(sample))
        names(tpm_data)[2] <- sample_id
        
        all_rna_raw <- left_join(all_rna_raw, tpm_data, by=c("gene_name"="gene_name"))
        
      }
      
    } else {
      cat("rna not found for ",sample_id, "\n")
      
    }
    
  }
  
  
  
  write.table( joined_file.filt,  file = file.path(secdir,"/results/rna.htseq.all.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  write.table( all_rna_raw.filt,  file = file.path(secdir,"/results/rna.tpm.all.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

} else {

  joined_file <- fread(file.path(secdir,"/results/rna.htseq.all.txt"))
  all_rna_raw <- fread(file.path(secdir,"/results/rna.tpm.all.txt"))

}

#### clean and test ####
## QC step, removing low mapping read samples
median_expression <- data.frame('med_ex' = apply(joined_file[,-1], 2, median, na.rm = TRUE))
ggplot(median_expression, aes(x=med_ex)) + geom_histogram() + theme_bw()

zero_median_samples <- names(which(apply(joined_file[,-1], 2, median, na.rm = TRUE) < 1))
print(zero_median_samples)

joined_file <- joined_file %>% dplyr::select( -all_of(zero_median_samples))
all_rna_raw <- all_rna_raw %>% dplyr::select( -all_of(zero_median_samples))

key_tumour_features <- fread(file.path(basedir,'/source_data/BTC_clinical.tumours.tsv'))
ffpe_tumours <- na.omit(key_tumour_features$Sample[key_tumour_features$`Sample Source(FFPE, FF)` == 'FFPE' ])

ffpe_tumours <- intersect(ffpe_tumours, names(all_rna_raw)[-1])

all_rna.ff <- dplyr::select(all_rna_raw, -all_of(ffpe_tumours))
rna.htseq.ff <- dplyr::select(joined_file, -all_of(ffpe_tumours))

#### manuscript cohort ####

overwrite.mn.cohort = T
if(overwrite.mn.cohort){
  
  good_sample_list <- fread(file.path(basedir,'/source_data/sample.list.txt'), header=F)
  
  tpm.good <- dplyr::select(all_rna.ff, all_of(c('gene_name',intersect(good_sample_list$V2,names(all_rna.ff)))))
  htseq.good <- dplyr::select(rna.htseq.ff, all_of(c('gene_name',intersect(good_sample_list$V2,names(rna.htseq.ff)))))
  
  
  
  write.table( tpm.good,  file = file.path(basedir,"/results/rna.tpm.ff.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  write.table( htseq.good,  file = file.path(secdir,"/rna.read.ff.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
}


```
