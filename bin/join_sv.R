#! /usr/bin/env Rscript

library(dplyr)

'%ni%' <- function(x,y)!('%in%'(x,y))

summarize_caller_pulls <- function(this.sv, this.caller){
  pull_this <- c('call'=NA,'vaf'=NA)
  if(this.caller %in% this.sv$CALLER){
    pull_this['call'] = this.sv$ALT[this.sv$CALLER == this.caller]
    if(length(pull_this['call']) > 1){pull_this['call'] = 'MULTIPLE'}
    pull_this['vaf'] = round(mean(this.sv$vaf[this.sv$CALLER == this.caller]), 2)
  } 
  return(pull_this)
}

BUFFER_SIZE = 10

library(optparse)

option_list = list(
  make_option(c("-t", "--tumour_id"), type="character", default=NULL, help="tumour_id", metavar="character"),
  make_option(c("-d", "--donor"), type="character", default=NULL, help="donor", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, help="output file path", metavar="character")
  
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

output_directory = opt$output
cat('output directory set to', output_directory, '\n')
tumour_id = opt$tumour_id
donor = opt$donor

ROOT_DIR = '/.mounts/labs/PCSI/'

##test
#tumour_id = 'BTC_0147_Lv_P_526'
#donor = 'BTC0147'
#ROOT_DIR = '/Volumes/pcsi/'


hg38_refSeq_genes <- read.table(paste0(ROOT_DIR,'/users/fbeaudry/hg38_refSeq_genes.txt'), header = T)
manta.path = paste0(ROOT_DIR,'pipeline/hg38_production/',donor,'/',tumour_id,'/wgs/bwa/0.7.17/manta/1.6.0/analysis/results/variants/somaticSV.vcf.gz')

delly.path = paste0(ROOT_DIR,'pipeline/hg38_production/',donor,'/',tumour_id,'/wgs/bwa/0.7.17/delly2/1.0.3/filter/',tumour_id,'_delly2_filter.vcf')

if( (file.exists(delly.path) | file.exists(paste0(delly.path, '.gz'))) ){
  cat('\ndelly2 1.0.3 exists\n')
} else {
  cat('\nsomeone messed with the file system\n')
  delly.path = paste0(ROOT_DIR,'pipeline/hg38_production/',donor,'/',tumour_id,'/wgs/bwa/0.7.17/delly2/1.0.3.q10/filter/',tumour_id,'_delly2_filter.vcf')
  
  if( (file.exists(delly.path) | file.exists(paste0(delly.path, '.gz'))) ){
    cat('\ndelly2 1.0.3.q10 exists\n')
  } else {
    cat('\nsomeone really messed with the file system\n')
    delly.path = paste0(ROOT_DIR,'pipeline/hg38_production/',donor,'/',tumour_id,'/wgs/bwa/0.7.17/delly2/1.0.3.q10s15z5/filter/',tumour_id,'_delly2_filter.vcf')
    
  }
  
}

if(  (file.exists(delly.path) | file.exists(paste0(delly.path, '.gz')))   ){
  #### delly ####
  cat('delly file found, filtering\n')
  INCLUDE_DELLY = TRUE
  
  ##delly
  
  if(file.exists(delly.path)){
    delly.data <- read.table(delly.path)
  } else if(file.exists(paste0(delly.path,'.gz'))){
    delly.data <- read.table(gzfile(paste0(delly.path,'.gz')))
  } 
  
  
  names(delly.data) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','tumour','normal')
  
  first=T
  for(this.index in c(1:nrow(delly.data))){
    sv_type = ''
    sv_length = NA
    this_row <- delly.data[this.index,]
    second_break = NULL
    
    this_info_split <- strsplit(this_row$INFO, ';')[[1]]
    for(this_info_index in c(1:length(this_info_split))){
      this_info_with_name <- strsplit(this_info_split[this_info_index], '=')[[1]]
      if(length(this_info_with_name) > 1){
        if(this_info_with_name[1] == 'SVTYPE'){
          sv_type = this_info_with_name[2]
        } else if(this_info_with_name[1] == 'END'){
          sv_length = abs(as.numeric(this_info_with_name[2]) - this_row$POS)
          second_break = as.numeric(this_info_with_name[2])
        }
        
      }
    }

    
    this_genotype_names <- strsplit(this_row$FORMAT, ':')[[1]]
    this_genotype_split <- strsplit(this_row$tumour, ':')[[1]]
    
    ##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference pairs">
    ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant pairs">
    
    if(is.null(second_break)){
      
    delly.processed <- cbind.data.frame(
      'CHROM' = this_row$CHROM,
      'POS' = this_row$POS,
      'REF' = this_row$REF,
      'ALT' = this_row$ALT,
      'SV_TYPE' = sv_type,
      'LENGTH' = sv_length,
      'CALLER' = 'DELLY',
      'VARIANT_SUPPORT' = as.numeric(this_genotype_split[this_genotype_names == 'DV']) ,
      'TOTAL_SUPPORT' = 
       ( as.numeric(this_genotype_split[this_genotype_names == 'DR']) + 
           as.numeric(this_genotype_split[this_genotype_names == 'DV']))
    )
    
    if(first==T){sv.calls <- delly.processed; first=F}else{
      sv.calls <- rbind.data.frame(sv.calls, delly.processed)
    }
    
    } else {

      if ( grepl(">", this_row$ALT) ) {
        this.alt <- paste0("]",this_row$CHROM, ":",second_break,"]")
      } else {
        this.alt <- this_row$ALT
      }
      
        delly.processed <- cbind.data.frame(
          'CHROM' = this_row$CHROM,
          'POS' = this_row$POS,
          'REF' = this_row$REF,
          'ALT' = this.alt,
          'SV_TYPE' = sv_type,
          'LENGTH' = sv_length,
          'CALLER' = 'DELLY',
          'VARIANT_SUPPORT' = as.numeric(this_genotype_split[this_genotype_names == 'DV']) ,
          'TOTAL_SUPPORT' = 
            ( as.numeric(this_genotype_split[this_genotype_names == 'DR']) + 
                as.numeric(this_genotype_split[this_genotype_names == 'DV']))
        )
        
        if(first==T){sv.calls <- delly.processed; first=F}else{
          sv.calls <- rbind.data.frame(sv.calls, delly.processed)
        }
        
        if ( grepl(">", this_row$ALT) ) {
          this.alt <- paste0("[",this_row$CHROM, ":",this_row$POS,"[")
        } else {
          this.alt <- this_row$ALT
        }
    
      delly.processed <- cbind.data.frame(
        'CHROM' = this_row$CHROM,
        'POS' = second_break,
        'REF' = this_row$REF,
        'ALT' = this.alt,
        'SV_TYPE' = sv_type,
        'LENGTH' = sv_length,
        'CALLER' = 'DELLY',
        'VARIANT_SUPPORT' = as.numeric(this_genotype_split[this_genotype_names == 'DV']) ,
        'TOTAL_SUPPORT' = 
          ( as.numeric(this_genotype_split[this_genotype_names == 'DR']) + 
              as.numeric(this_genotype_split[this_genotype_names == 'DV']))
      )
      sv.calls <- rbind.data.frame(sv.calls, delly.processed)
    }
    
  }
}

if( file.exists(manta.path)  ){
  #### manta ####
  INCLUDE_MANTA = TRUE
  cat('manta file found, filtering\n')
  
  manta.data <- read.table(gzfile(manta.path))
  names(manta.data) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','normal','tumour')
  
  manta.data <- manta.data %>% filter(FILTER == 'PASS')
  
  for(this.index in c(1:nrow(manta.data))){
    
    sv_type = ''
    sv_length = NA
    this_row <- manta.data[this.index,]
    second_break = NULL
    
    this_info_split <- strsplit(this_row$INFO, ';')[[1]]
    
    for(this_info_index in c(1:length(this_info_split))){
      this_info_with_name <- strsplit(this_info_split[this_info_index], '=')[[1]]
      if(length(this_info_with_name) > 1){
        
        if(this_info_with_name[1] == 'SVTYPE'){
          sv_type = this_info_with_name[2]
        }  else if(this_info_with_name[1] == 'END'){
          sv_length = abs(as.numeric(this_info_with_name[2]) - this_row$POS)
          second_break = as.numeric(this_info_with_name[2])
        }
        
      }
    }
    
    ##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Spanning paired-read support for the ref and alt alleles in the order listed">
    ##FORMAT=<ID=SR,Number=.,Type=Integer,Description="Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999">
    this_genotype_names <- strsplit(this_row$FORMAT, ':')[[1]]
    this_genotype_split <- strsplit(this_row$tumour, ':')[[1]]
    if("SR" %in% this_genotype_names){
      ref_alt <- strsplit(this_genotype_split[this_genotype_names == 'SR'], ',')[[1]]
    } else {
      ref_alt <- strsplit(this_genotype_split, ',')[[1]]
    }
    
    if(is.null(second_break)){
      sv.calls <- rbind.data.frame(sv.calls, 
         cbind.data.frame(
           'CHROM' = this_row$CHROM,
           'POS' = this_row$POS,
           'REF' = this_row$REF,
           'ALT' = this_row$ALT,
           'SV_TYPE' = sv_type,
           'LENGTH' = sv_length,
           'CALLER' = 'MANTA',
           'VARIANT_SUPPORT' = as.numeric(ref_alt[2]) ,
           'TOTAL_SUPPORT' = 
             ( as.numeric(ref_alt[1]) + 
                 as.numeric(ref_alt[2]))
         )
      )
    }  else {
        
        if ( grepl(">", this_row$ALT) ) {
          this.alt <- paste0("]",this_row$CHROM, ":",second_break,"]")
        } else {
          this.alt <- this_row$ALT
        }
        
        sv.calls <- rbind.data.frame(sv.calls, 
             cbind.data.frame(
               'CHROM' = this_row$CHROM,
               'POS' = this_row$POS,
               'REF' = this_row$REF,
               'ALT' = this.alt,
               'SV_TYPE' = sv_type,
               'LENGTH' = sv_length,
               'CALLER' = 'MANTA',
               'VARIANT_SUPPORT' = as.numeric(ref_alt[2]) ,
               'TOTAL_SUPPORT' = 
                 ( as.numeric(ref_alt[1]) + 
                     as.numeric(ref_alt[2]))
             )
        )
        
        if ( grepl(">", this_row$ALT) ) {
          this.alt <- paste0("[",this_row$CHROM, ":",this_row$POS,"[")
        } else {
          this.alt <- this_row$ALT
        }
        
        
      sv.calls <- rbind.data.frame(sv.calls, 
         cbind.data.frame(
           'CHROM' = this_row$CHROM,
           'POS' = second_break,
           'REF' = this_row$REF,
           'ALT' = this.alt ,
           'SV_TYPE' = sv_type,
           'LENGTH' = sv_length,
           'CALLER' = 'MANTA',
           'VARIANT_SUPPORT' = as.numeric(ref_alt[2]) ,
           'TOTAL_SUPPORT' = 
             ( as.numeric(ref_alt[1]) + 
                 as.numeric(ref_alt[2]))
         )
      )
    }
    
  }
}

svaba.path = paste0(ROOT_DIR,'pipeline/hg38_production/',donor,'/',tumour_id,'/wgs/bwa/0.7.17/SVaBA/v134/withdbSNP/',tumour_id,'.svaba.somatic.sv.addsvtype.vcf')
new.svaba.path = paste0(ROOT_DIR,'pipeline/hg38_production/',donor,'/',tumour_id,'/wgs/bwa/0.7.17/SVaBA/1.2.0/withdbSNP/',tumour_id,'.svaba.somatic.sv.addsvtype.vcf.gz')

if(  (file.exists(svaba.path) | file.exists(paste0(svaba.path, '.gz')) |  file.exists(new.svaba.path) )    ){
  #### svaba ####
  INCLUDE_SVABA = TRUE
  cat('svaba file found, filtering\n')
  
  if(file.exists(svaba.path)){
    svaba.data <- read.table(svaba.path)
  } else if(file.exists(new.svaba.path)){
    svaba.data <- read.table(new.svaba.path)
  } else {
    svaba.data <- read.table(gzfile(paste0(svaba.path,'.gz')))
  }
    
  names(svaba.data) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','normal','tumour')
  
  for(this.index in c(1:nrow(svaba.data))){
    sv_type = ''
    sv_length = NA
    this_row <- svaba.data[this.index,]
    
    this_info_split <- strsplit(this_row$INFO, ';')[[1]]
    for(this_info_index in c(1:length(this_info_split))){
      this_info_with_name <- strsplit(this_info_split[this_info_index], '=')[[1]]
      if(length(this_info_with_name) > 1){
        if(this_info_with_name[1] == 'SVTYPE_ID'){
          sv_type = this_info_with_name[2]
        }else if(this_info_with_name[1] == 'SPAN'){
          sv_length = as.numeric(this_info_with_name[2])
        }
        
      }
    }
    
    this_genotype_names <- strsplit(this_row$FORMAT, ':')[[1]]
    this_genotype_split <- strsplit(this_row$tumour, ':')[[1]]
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth of coverage: Number of reads covering site.">
    ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allele depth: Number of reads supporting the variant">
    

    
    sv.calls <- rbind.data.frame(sv.calls, 
                                 cbind.data.frame(
                                   'CHROM' = this_row$CHROM,
                                   'POS' = this_row$POS,
                                   'REF' = this_row$REF,
                                   'ALT' = this_row$ALT,
                                   'SV_TYPE' = sv_type,
                                   'LENGTH' = sv_length,
                                   'CALLER' = 'SVABA',
                                   'VARIANT_SUPPORT' = as.numeric(this_genotype_split[this_genotype_names == 'AD']) ,
                                   'TOTAL_SUPPORT' =  as.numeric(this_genotype_split[this_genotype_names == 'DP'])
                                 )
    )
    
    
  }
}

sv.sorted <- sv.calls %>% arrange(CHROM, POS )

sv.sorted$vaf <- sv.sorted$VARIANT_SUPPORT / sv.sorted$TOTAL_SUPPORT
sv.sorted$vaf[ sv.sorted$vaf > 1] <- 1
sv.sorted$vaf[ sv.sorted$TOTAL_SUPPORT == 0] <- 0

cat('printing raw concatenated data to', paste0(output_directory, tumour_id,'.sv.txt'), '\n')
write.table(
  sv.sorted,
  file = paste0(output_directory, tumour_id,'.sv.txt'),
  row.names = FALSE, quote = FALSE, sep = "\t"
)


#### start filtering ####
cat('starting combining\n')

#sv.sorted <- read.table(paste0(output_directory,tumour_id,'.sv.txt'), header = T)

sv.filtered <- sv.sorted
## regroup variants into single rows

row.index=1
while(row.index <= nrow(sv.filtered)){
  
  this.start <- sv.filtered$POS[row.index]
  this.chr <- sv.filtered$CHROM[row.index]
  
  this.end <- this.start + BUFFER_SIZE
  this.sv <- sv.filtered[sv.filtered$POS >= this.start & sv.filtered$POS <= this.end & sv.filtered$CHROM == this.chr,]

  #gonna assume this extension of the window needs only happen once
  this.end <- max(this.sv$POS) + BUFFER_SIZE
  this.sv <- sv.filtered[sv.filtered$POS >= this.start & sv.filtered$POS <= this.end & sv.filtered$CHROM == this.chr,]
  
  this.sv_type <- this.sv$SV_TYPE

  this.sv_type <- paste(unique(this.sv_type[this.sv_type != "BND"]), sep = ' ,')
  if(length(this.sv_type) == 0){this.sv_type = "BND"}

  this.support = 0
  if(!is.na(summarize_caller_pulls(this.sv, 'DELLY')['call'])){this.support = this.support + 1}
  if(!is.na(summarize_caller_pulls(this.sv, 'MANTA')['call'])){this.support = this.support + 1}
  if(!is.na(summarize_caller_pulls(this.sv, 'SVABA')['call'])){this.support = this.support + 1}
  
  
  this.sv_summary <- 
  cbind.data.frame(
    'CHROM' = unique(this.sv$CHROM),
    'POS' = as.numeric(names(sort(table(this.sv$POS),decreasing=TRUE))[1]),
    'REF' = names(sort(table(this.sv$REF),decreasing=TRUE))[1],
    'SV' = this.sv_type,
    #svaba has weird VAFs > 1
    'VAF' = round(mean(this.sv$vaf[this.sv$CALLER %ni% c('SVABA')]), 3),
    'LENGTH' = round(mean(this.sv$LENGTH), 0),
    'COV' = round(mean(this.sv$TOTAL_SUPPORT[this.sv$CALLER %ni% c('SVABA')]), 2),
    'SUPPORT' = this.support,
    'DELLY_CALL' = summarize_caller_pulls(this.sv, 'DELLY')['call'],
    'DELLY_VAF' = summarize_caller_pulls(this.sv, 'DELLY')['vaf'],
    'MANTA_CALL' = summarize_caller_pulls(this.sv, 'MANTA')['call'],
    'MANTA_VAF' = summarize_caller_pulls(this.sv, 'MANTA')['vaf'],
    'SVABA_CALL' = summarize_caller_pulls(this.sv, 'SVABA')['call'],
    'SVABA_VAF' = summarize_caller_pulls(this.sv, 'SVABA')['vaf']
  )
  
  if(row.index==1){sv.summary<- this.sv_summary}else{sv.summary<- rbind.data.frame(sv.summary, this.sv_summary)}
    
    
  row.index = row.index + nrow(this.sv)
  
}

sv.summary$SV_CLASS[sv.summary$SV %in% c("DEL", "DUP")] <- 'UNBALANCED'
sv.summary$SV_CLASS[sv.summary$SV %in% c("INV", "TRA")] <- 'BALANCED'

sv.summary$GENE <- NA
for(row.index in c(1:nrow(sv.summary))){
  this_gene <- unique(hg38_refSeq_genes$gene[hg38_refSeq_genes$chrom == sv.summary$CHROM[row.index] &
                                               hg38_refSeq_genes$start <= sv.summary$POS[row.index] &
                                               hg38_refSeq_genes$end >= sv.summary$POS[row.index]])
  if(length(this_gene) > 0){
    sv.summary$GENE[row.index] <- this_gene
  }
}


write.table(
  sv.summary,
  file = paste0(output_directory, tumour_id, '.sv.filtered.txt'),
  row.names = FALSE, quote = FALSE, sep = "\t"
)
sv.summary.saved <- sv.summary

#### ADD ORTHOGONAL ####

#sv.summary <- read.table(paste0(output_directory,tumour_id,'.sv.filtered.txt'), header = T)

celluloid_solution = 'solution'
#celluloid_solution = 'solution1'
celluloid.path = paste0(ROOT_DIR,'/pipeline/hg38_production/',donor,'/',tumour_id,'/wgs//bwa/0.7.17/celluloidXY/v0.11.7/',celluloid_solution,'/segments_',tumour_id,'.txt.sorted.bed.gz')

if(file.exists(celluloid.path) ){
  
  #### celluloid ####
  INCLUDE_CELLULOID = TRUE
  
  celluloid.data <- read.table(gzfile(celluloid.path))
  names(celluloid.data) <- c('chrom','start.pos','end.pos','size','arm','sampleID','mean','meanmap','mask','meanar','p','imean','labels')
  
  celluloid.data <- celluloid.data %>% filter(mask == 'FALSE')
  
  unbalanced.sv.summary <- sv.summary %>% filter( CHROM %in% unique(celluloid.data$chrom) & SUPPORT >1) #SV_CLASS == 'UNBALANCED' &
  unbalanced.sv.summary$CELLULOID_POS <- NA
  unbalanced.sv.summary$CELLULOID_CALL <- NA
  
  CELLULOID_BUFFER = 100000

  for(this.chrom in unique(celluloid.data$chrom)){
    
    this.celluloid.data <- celluloid.data[celluloid.data$chrom == this.chrom,]
    this.sv.summary <- unbalanced.sv.summary[unbalanced.sv.summary$CHROM == this.chrom ,]
    
    these_cnv_pos = c()
    these_sv_pos = c()
    
    if(nrow(this.sv.summary) > 0){
      
      for(this.cnv.pos in this.celluloid.data$start.pos ){
        
        for(this.sv.pos in this.sv.summary$POS ){
        
          if( ( (this.sv.pos + CELLULOID_BUFFER) > this.cnv.pos &
             (this.sv.pos - CELLULOID_BUFFER) <= this.cnv.pos ) ) {
            
                if(this.cnv.pos %in% unbalanced.sv.summary$CELLULOID_POS){
                  unbalanced.sv.summary <- unbalanced.sv.summary %>% filter(!(CELLULOID_POS == this.cnv.pos & CHROM == this.chrom & is.na(SV)))
                }
            
                unbalanced.sv.summary$CELLULOID_POS[unbalanced.sv.summary$CHROM == this.chrom & unbalanced.sv.summary$POS == this.sv.pos] <- this.cnv.pos
                unbalanced.sv.summary$CELLULOID_CALL[unbalanced.sv.summary$CHROM == this.chrom & unbalanced.sv.summary$POS == this.sv.pos] <-  round(this.celluloid.data$imean[this.celluloid.data$start.pos == this.cnv.pos], 2)
            
                these_cnv_pos = c(these_cnv_pos, this.cnv.pos)
                these_sv_pos = c(these_sv_pos, this.sv.pos)
                
          } else if(this.cnv.pos %ni% these_cnv_pos & this.sv.pos %ni% these_sv_pos){
            
            these_cnv_pos = c(these_cnv_pos, this.cnv.pos)
            these_sv_pos = c(these_sv_pos, this.sv.pos)

              this_gene <- unique(hg38_refSeq_genes$gene[hg38_refSeq_genes$chrom == this.chrom &
                                                           hg38_refSeq_genes$start <= this.cnv.pos &
                                                           hg38_refSeq_genes$end >= this.celluloid.data$end.pos[this.celluloid.data$start.pos == this.cnv.pos]])
              if(length(this_gene) == 0){
                 this_gene <- NA
              }
              if(length(this_gene) > 1){
                this_gene <- paste(unlist(this_gene))
              }
              
            #cat(this.chrom, this.cnv.pos, this.sv.pos, this_gene, '\n')  
            unbalanced.sv.summary <- rbind.data.frame(
              unbalanced.sv.summary,
              cbind.data.frame(
                'CHROM' = this.chrom,
                'POS' = this.cnv.pos,
                'REF' = NA,
                'SV' = NA,
                #svaba has weird VAFs > 1
                'VAF' = 1,
                'LENGTH' = this.celluloid.data$size[this.celluloid.data$start.pos == this.cnv.pos],
                'COV' = NA,
                'SUPPORT' = 1,
                'DELLY_CALL' = NA,
                'DELLY_VAF' = NA,
                'MANTA_CALL' = NA,
                'MANTA_VAF' = NA,
                'SVABA_CALL' = NA,
                'SVABA_VAF' = NA,
                'SV_CLASS' = 'UNBALANCED',
                'GENE' =this_gene,
                'CELLULOID_POS' = this.cnv.pos,
                'CELLULOID_CALL' = round(this.celluloid.data$imean[this.celluloid.data$start.pos == this.cnv.pos], 2)
                
              )
            )
            
            
            
            }
         }
        }
      }
    }  
  }
  
sv.summary.sorted <- unbalanced.sv.summary[order(unbalanced.sv.summary$CHROM, unbalanced.sv.summary$POS ),]


sv.cnv <- sv.summary.sorted

this.sv.count <- nrow(sv.cnv) 
this.expected.sv <- this.sv.count / 3000000000

this.bin.size <- 3000000000 /  this.sv.count 
Expected=1

sv.cnv$OE2E <- 0
for(row.index in c(1:nrow(sv.cnv))){
  
  bin.start = sv.cnv$POS[row.index] - this.bin.size
  bin.end = sv.cnv$POS[row.index] + this.bin.size
  
  sv.bin <- sv.cnv[sv.cnv$POS >bin.start & sv.cnv$POS < bin.end & sv.cnv$CHROM == sv.cnv$CHROM[row.index],]
  sv.observed <- nrow(sv.bin)
  OE2E <- (sv.observed-Expected)^2 / Expected
  sv.cnv$OE2E[row.index] <- OE2E
  
  ##TODO: if in bin where OE2E > 100, type = 'CLUSTER'
  if(OE2E > 100){
    sv.cnv$SV[sv.cnv$POS >bin.start & sv.cnv$POS < bin.end & sv.cnv$CHROM == sv.cnv$CHROM[row.index]] <- 'CLU'
    
  }
  
}
sv.cnv$tumour <- tumour_id

write.table(
  sv.cnv,
  file = paste0(output_directory, tumour_id, '.sv.cnv.txt'),
  row.names = FALSE, quote = FALSE, sep = "\t"
)


#### VAF for BALANCED SVs ####

celluloid.param.path = paste0(ROOT_DIR,'/pipeline/hg38_production/',donor,'/',tumour_id,'/wgs//bwa/0.7.17/celluloidXY/v0.11.7/',celluloid_solution,'/parameters_',tumour_id,'.txt')

#sv.summary <- read.table(paste0(output_directory, tumour_id, '.sv.filtered.txt'), header = T)
sv.summary <- sv.summary.saved
sv.balanced <- sv.summary %>% filter(SV_CLASS == 'BALANCED' & SUPPORT > 1)

sv.balanced$major_cn <- NA
sv.balanced$minor_cn <- NA
for(this.sv.index in c(1:nrow(sv.balanced)) ){
  this_sv <- sv.balanced[this.sv.index,]
  cn.alleles <- celluloid.data$labels[celluloid.data$chrom == this_sv$CHROM & 
                          celluloid.data$start.pos <= this_sv$POS &
                   celluloid.data$end.pos > this_sv$POS ]
  if(length(cn.alleles) > 0 ){
    split.alleles <- strsplit(x = as.character(cn.alleles), split = '\\.')
    sv.balanced$major_cn[this.sv.index] <- split.alleles[[1]][2]
    sv.balanced$minor_cn[this.sv.index] <- split.alleles[[1]][1]
  }
}

sv.balanced$mutation_id <- paste(sv.balanced$CHROM, sv.balanced$POS, sv.balanced$SV, sep=':')
sv.balanced$sample_id <- 'sample1'
sv.balanced$alt_counts <- round(sv.balanced$VAF * sv.balanced$COV, 0)
sv.balanced$ref_counts <- round(sv.balanced$COV, 0) - sv.balanced$alt_counts
sv.balanced$normal_cn <- 2

celluloid.param <- read.table(celluloid.param.path, header = T)
sv.balanced$tumour_content <- round(celluloid.param$T1, 2)

##mutation_id	sample_id	ref_counts	alt_counts	normal_cn	major_cn	minor_cn	tumour_content
##chr1:3918869:AT	sample1	48	23	2	1	0	0.67

sv.pyclone <- sv.balanced %>% 
  dplyr::select(mutation_id,sample_id,ref_counts,alt_counts,normal_cn, major_cn, minor_cn, tumour_content) %>% 
  filter(!is.na(major_cn) & !is.na(minor_cn))

write.table(
  sv.pyclone,
  file = paste0(output_directory, tumour_id, '.sv.depths.txt'),
  row.names = FALSE, quote = FALSE, sep = "\t"
)

#mutation_id	rmsk	simplerepeat	superdup	ANNOVAR_EXONIC	ANNOVAR
sv.balanced.info <- sv.balanced %>% 
  dplyr::select(mutation_id,GENE)

write.table(
  sv.balanced.info,
  file = paste0(output_directory, tumour_id, '.sv.info.txt'),
  row.names = FALSE, quote = FALSE, sep = "\t"
)

#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	TUMOR
#chr1	2054566	chr1:2054566:G	GT	G	.	PASS	t_alt_count=0;t_ref_count=0	GT:AD:DP:FT	0/1:0,0:0:PASS

sv.balanced.vcf.raw <- sv.balanced
sv.balanced.vcf.raw$QUAL <- '.'
sv.balanced.vcf.raw$FILTER <- 'PASS'
sv.balanced.vcf.raw$FORMAT <- 'GT:AD:DP:FT'
sv.balanced.vcf.raw$INFO <- paste0('t_alt_count=',sv.balanced.vcf.raw$alt_counts, ';t_ref_count=',sv.balanced.vcf.raw$ref_counts)
sv.balanced.vcf.raw$TUMOR <- paste0('0/1:',sv.balanced.vcf.raw$ref_counts,',',sv.balanced.vcf.raw$alt_counts,':',round(sv.balanced.vcf.raw$COV,0),':','PASS')

sv.balanced.vcf <- sv.balanced.vcf.raw %>% dplyr::select(CHROM, POS, mutation_id, REF, SV, QUAL, FILTER, INFO, FORMAT, TUMOR)

write.table(
  sv.balanced.vcf,
  file = paste0(output_directory, tumour_id, '.sv.mutationTimeR.headless.vcf'),
  row.names = FALSE, quote = FALSE, sep = "\t", col.names = F
)
  







