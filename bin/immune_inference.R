```{r}

suppressPackageStartupMessages({

library(data.table)
library(dplyr)
library(immunedeconv)
  
})

basedir = '~/Documents/scripts/github/btc-ccs/'
secdir = '~/Documents/scripts/github/evolution/'

set_cibersort_binary(paste0(secdir, '/code/src/CIBERSORT.R'))

# Cibersort loadings matrix file
set_cibersort_mat(paste0(basedir, '/source_data/LM22.txt')) 

all_rna <- fread(file.path(basedir,"/results/rna.tpm.ff.txt"), header = T)

TPMs_full <- as.data.frame(all_rna[,-1])
rownames(TPMs_full) <- all_rna$gene_name

results_dataframe = immunedeconv::deconvolute(TPMs_full, "cibersort_abs")

write.table( results_dataframe,  
             file = file.path(basedir,"/results/LBR.cibersort.txt"), 
             row.names = F, quote = FALSE, sep = "\t")

```
