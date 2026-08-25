# Integrated whole-genome and transcriptome sequencing reveals divergent evolutionary processes across biliary tract cancer subtypes
Discover and compare the consensus cancer subtypes of biliary tract cancer for manuscript under review (and on [bioRxiv](https://www.biorxiv.org/content/10.64898/2025.12.12.693962v3)).

## Quick Classification
For quick classification of a bulk tumor specimen into expression subtypes, run:

`Rscript call_btc_eCCS.R -i data/tpm.txt -o data/eCCS.txt -m data/LBR.tps.classifier.rds` 

where input.txt is a table of RNA expression data (trained with TPM but anecdotally robust to other quantification metrics) with the first column as gene identifiers. 

For gCCS classification, you'll need CNTools (`BiocManager::install("CNTools")`). Then run:

`Rscript call_btc_gCCS.R -s data/example.seg -v data/mean_vaf.txt -o data/gCCS.txt -m data/cnv.model.rds -b data/genebed.hg38.txt -c data/cytoband.txt`

## Results walkthrough 
To see how we generated our results, check out the full analysis of the data on our [r-markdown generated site](https://pancurx-oicr.github.io/btc-ccs/) showing step by step walkthroughs of all analyses using the data and code in this repository. 

## Start from scratch 
If you are looking instead to start from scratch, raw .fastq files can be found on EGA under study ID [EGAS50000000972](https://ega-archive.org/studies/EGAS50000000972). Variants calls are generated according to the pipeline of [Chan-Seng-Yue 2020](https://www.nature.com/articles/s41588-019-0566-9). HPC analyses are launched with `0_evolution.collations.sh`, for which shell scripts are found in `bin/`. Data can be accessed with a Data Access Agreement.

## License and Usage Restrictions

### Code
Copyright © 2026 Felix Beaudry.  
All rights reserved.

This repository and all source code contained herein are **unpublished, proprietary work** and are provided for **viewing purposes only**.  

No part of this code may be copied, modified, distributed, or used, in whole or in part, without prior written permission from the author.
