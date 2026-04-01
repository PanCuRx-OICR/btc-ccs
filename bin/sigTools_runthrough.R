#! /usr/bin/env Rscript

library(optparse)
library(tidyr)

option_list = list(
  make_option(c("-t", "--tumour_id"), type="character", default=NULL, help="tumour_id", metavar="character"),
  make_option(c("-S", "--snvFile"), type="character", default=NULL, help="SNV vcf file", metavar="character"),
  make_option(c("-I", "--indelFile"), type="character", default=NULL, help="indel vcf file", metavar="character"),
  make_option(c("-V", "--SVFile"), type="character", default=NULL, help="Structural variant file", metavar="character"),
  make_option(c("-L", "--LOHFile"), type="character", default=NULL, help="LOH file", metavar="character"),
  make_option(c("-g", "--genomeVersion"), type="character", default="hg38", help="genome version", metavar="character"),
  make_option(c("-o", "--outpath"), type="character", default="hg38", help="genome version", metavar="character"),
  make_option(c("-b", "--basedir"), type="character", default=NULL, help="location of call_hrdetect", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)


genomeVersion       <-  opt$genomeVersion
sbs_location     <-  opt$snvFile
indel_vcf_location  <-  opt$indelFile
SV_vcf_location     <-  opt$SVFile
LOH_seg_file        <-  opt$LOHFile
base_dir            <-  opt$basedir
sample_name         <-  opt$tumour_id
outpath             <-  opt$outpath

source(paste0(base_dir, "/call_hrdetect.R"))

if(file.exists(indel_vcf_location) & file.exists(LOH_seg_file) & file.exists(SV_vcf_location) & file.exists(sbs_location) ){
  
  input_matrix <- matrix(NA, nrow = 1, ncol = 6,
                         dimnames = list(
                           sample_name,
                           c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
                         )
  )
    
#### pre-process SEGMENT data ####
  
  cat('pre-processing LOH data\n')
  seg.data <- read.table(LOH_seg_file, sep="\t", header=FALSE)
  
  LOH_ls    <- summarize_LOH(seg.data)
  input_matrix[sample_name,"hrd"]  <- unlist(LOH_ls[1])
  
#### pre-process In/Del data ####
  cat('pre-processing in/del data\n')
  
  expected_chroms <- paste0("chr", c(seq(1:22), "X", "Y"))
  genomeSeq <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  gr <- GenomicRanges::GRanges(GenomeInfoDb::seqinfo(genomeSeq))
  
  indel_vcf <- read.table(indel_vcf_location, sep = "\t", header = TRUE)
  vcf_seqnames <- Rsamtools::headerTabix(indel_vcf_location)$seqnames
  gr <- GenomeInfoDb::keepSeqlevels(gr, intersect(vcf_seqnames, 
                                                  expected_chroms), pruning.mode = "coarse")
  GenomeInfoDb::seqlevels(gr) <- sub("chr", "", GenomeInfoDb::seqlevels(gr))
  indel.data <- VariantAnnotation::readVcf(file = indel_vcf_location, 
                                           genome = "hg38" )
  
  indel.data <- VariantAnnotation::expand(indel.data)
  
  indel.df <- prepare.indel.df(indel.data, genomeSeq, genome.v=genomeVersion, 
                               expected_chroms)
  

  indels_classified <- mh(indel.df)
  indel_count_proportion <- indelsToCountAndProportion(indels_classified, 
                                                     sample_name)
  

  indel_ls  <- summarize_indels(indel_count_proportion)
  input_matrix[sample_name,"del.mh.prop"] <- indel_ls$sim_proportions["del.mh.prop"]

#### pre-process Structural Variant data ####
  cat('pre-processing SV data\n')
  
  SV_vcf        <- read.table(SV_vcf_location, sep = "\t", header = TRUE)
  
  SV_vcf$length_bin <- 'tiny'
  SV_vcf$length_bin[SV_vcf$LENGTH > 1000 & SV_vcf$LENGTH <= 10000] <- '1-10kb'
  SV_vcf$length_bin[SV_vcf$LENGTH > 10000 & SV_vcf$LENGTH <= 100000] <- '10-100kb'
  SV_vcf$length_bin[SV_vcf$LENGTH > 100000 & SV_vcf$LENGTH <= 1000000] <- '100kb-1Mb'
  SV_vcf$length_bin[SV_vcf$LENGTH > 1000000 & SV_vcf$LENGTH <= 10000000] <- '1-10Mb'
  SV_vcf$length_bin[SV_vcf$LENGTH > 10000000 ] <- '>10Mb'
  
  sv_tally <- SV_vcf %>% group_by(length_bin, SV) %>% tally()

  ## 'non-clustered_tds_1-10Kb' + 'non-clustered_tds_10-100Kb' + 'non-clustered_trans'
  input_matrix[sample_name,"SV3"]  <- sum(sv_tally$n[sv_tally$SV == 'DUP' & sv_tally$length_bin %in% c('1-10kb', '10-100kb')], na.rm = T)
  
  ## 'non-clustered_del_10-100Kb + 'non-clustered_del_1-10Kb' + 'non-clustered_trans'
  input_matrix[sample_name,"SV5"]  <- sum(sv_tally$n[sv_tally$SV == 'DEL' & sv_tally$length_bin %in% c('1-10kb', '10-100kb')], na.rm = T)
  
#### pre-process SNV data ####
  cat('pre-processing SNV data\n')
  
  snv_df    <- read.table(sbs_location, header = TRUE)

  input_matrix[sample_name,"SNV3"] <- snv_df$Signature.3
  input_matrix[sample_name,"SNV8"] <- snv_df$Signature.8
  
  
#### Call HRDetect ####
  
  HRDetect_res <- applyHRDetectDavies2017(input_matrix)
  
  HRDetect_res <- cbind.data.frame(HRDetect_res, input_matrix)
  HRDetect_res$tumour <- rownames(HRDetect_res)

  write.table(
    HRDetect_res,
    file = outpath ,
    row.names = FALSE, quote = FALSE, sep = "\t"
  )

}


  