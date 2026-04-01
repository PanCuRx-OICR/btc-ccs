```{r}

#devtools::install_github("mg14/mg14")
#devtools::install_github("gerstung-lab/MutationTimeR")
#remotes::install_github("caravagn/mobster")

library(PNWColors)
library(data.table)
library(ggplot2)
library(tidyr)
library(dplyr)
library(mobster)
library(VariantAnnotation)
library(MutationTimeR)

packageVersion('mobster')
packageVersion('MutationTimeR')

'%ni%' <- function(x,y)!('%in%'(x,y))

working_dir = '/Volumes/pcsi/users/fbeaudry/more.data/'

script_dir='~/Documents/scripts/github/evolution/'
source(paste0(script_dir,'/MutationTimeR/R/utils.R'))
source(paste0(script_dir,'/MutationTimeR/R/MutationTime.R'))

ROOT_DIR = '/Volumes/pcsi/'
sample_list <- read.table(paste0(ROOT_DIR,'/users/fbeaudry/sample.list.txt'), header=F)

sv.codes = c('TRA','INV')
clonal.times = c("clonal [NA]", "clonal [late]",  "clonal [early]")

#### MOBSTER ####

these.samples <- sample_list$V2
for(tumour_id in these.samples){
  
  cat('processing ',tumour_id,'\n')
  
  working_dir = paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id)
  
  if (!dir.exists(file.path(working_dir, 'mobster'))) {
    dir.create(file.path(working_dir, 'mobster'), recursive = T)
  }
  

  out_dir = paste0(working_dir, '/mobster/')
  
  depths.file = paste0(working_dir,'/pyclone/',tumour_id,'.all.depths.txt')
  
  ####MOBSTER: analysis ####
  if(file.exists(depths.file)){
    depths <- fread(depths.file)
    
    #the +1 ensures VAF != 1, which mobster won't allow
    #depths$VAF <- (depths$alt_counts ) / (depths$alt_counts + depths$ref_counts + 1)
    
    depths$VAF <- (depths$alt_counts ) / (depths$alt_counts + depths$ref_counts)
    depths <- depths %>% dplyr::filter(VAF < 1 & VAF > 0)
    
    fit = mobster_fit(
      depths,    
      auto_setup = "FAST"    
      
    )
    
    mobster_clusters <- Clusters(fit$best) 
    possible_clusters <- unique(mobster_clusters$cluster)
    mobster_assignment_prob <- mobster_clusters %>% dplyr::select(possible_clusters)
    mobster_clusters$mobster_assignment_prob <- apply(mobster_assignment_prob, 1, max)
    mobster_clusters$sample_id <- tumour_id
    
    mobster_clusters <- mobster_clusters %>% dplyr::select(sample_id,mutation_id,  cluster,mobster_assignment_prob)
    
    ####MOBSTER: output ####
    
    png(paste0(out_dir,'/',tumour_id,'.mobster.png') ,width = 1000, height = 400) #, bg = "transparent"
    print(
      plot(fit$best) + theme_bw(base_size = 22) + labs(title='') + theme(panel.grid = element_blank())
    )
    dev.off()
    
    
    write.table(
      mobster_clusters,
      file = paste0(out_dir,'/',tumour_id,'.all.mobster.txt'),
      row.names = FALSE, quote = FALSE, sep = "\t"
    )
    
    cat('found: ',possible_clusters,'\n')
    if('Tail' %in% possible_clusters ){
      write.table(
        evolutionary_parameters(fit),
        file = paste0(out_dir,'/',tumour_id,'.mobster.fit.txt'),
        row.names = FALSE, quote = FALSE, sep = "\t"
      )
    }
  } else {
    cat('some files missing for ',tumour_id,'\n')
    
  }
}

#### mutationTimeR ####

working_dir = '/Volumes/pcsi/users/fbeaudry/more.data/'

these_samples <- sample_list$V2 #[150:169]

for(tumour_id in these_samples){
  
  vcfFileIn.snv         = paste0(working_dir,'/',tumour_id,'/pyclone/',tumour_id,'.snv.mutationTimeR.vcf')
  pyclone.file          = paste0(working_dir,'/',tumour_id,'/pyclone/',tumour_id,'.all.pyclone.tsv')
  
  # requires files are header and celluloid links are respected, collated on cluster in evolution.collations.sh
  cnv_file              = paste0(working_dir,'/',tumour_id,'/mutationtimer/segments_',tumour_id,'.txt')
  celluloid_params_file = paste0(working_dir,'/',tumour_id,'/mutationtimer/parameters_',tumour_id,'.txt')
  vcfFileIn.sv          = paste0(working_dir,'/',tumour_id,'/mutationtimer/',tumour_id,'.sv.mutationTimeR.vcf')
  
  if(file.exists(pyclone.file) & (file.exists(cnv_file) | file.exists(paste0(cnv_file,'.gz'))) & 
     file.exists(celluloid_params_file) & file.exists(vcfFileIn.snv) & file.exists(vcfFileIn.sv)){
    
    cat('processing ',tumour_id,'\n')
    
    
    pyclone <- fread(pyclone.file)
    
    clusters <- pyclone %>% group_by(cluster_id, cellular_prevalence) %>% tally()
    names(clusters) <- c("cluster","proportion", "n_ssms")
    clusters <- as.data.frame(clusters)
    clusters <- clusters[order(clusters$proportion, decreasing = T),]
    #renumber clusters after reordering
    clusters$cluster <- c(0:(length(clusters$cluster)-1))
    
    
    ##TIMER: cnv ##
    
    
    if(file.exists(paste0(cnv_file,'.gz'))){
      cnvs = read.table(gzfile(paste0(cnv_file,'.gz')), header = F)
    } else {
      cnvs = read.table(cnv_file, header = F)
    }
    cnvs <- cnvs %>% dplyr::filter(V1 != 'chrom')
    
    celluloid_params = fread(celluloid_params_file)
    purity = celluloid_params$T1
    
    # Copy number segments, needs columns  major_cn, minor_cn and clonal_frequency of each segment
    
    cnvs <- separate(data = cnvs, col = V13, into = c("minor_cn", "major_cn"), sep = "\\.")
    
    #this filter removes homozygous deletions (mutationTimeR can handle, but pyclone removes)
    cnvs <- cnvs[!is.na(cnvs$major_cn) & !is.na(cnvs$minor_cn),]
    
    
    #celluloid segments overlap by one base
    cnvs$V2 <- as.numeric(cnvs$V2)
    cnvs$V3 <- as.numeric(cnvs$V3)
    gr <- GRanges(seqnames = cnvs$V1, ranges = IRanges(cnvs$V2, cnvs$V3-1), strand = '*') 
    
    #setting clonal_frequency=purity here assumes no subclonal CNVs
    cn <- GRanges(gr, major_cn=as.integer(cnvs$major_cn) , minor_cn=as.integer(cnvs$minor_cn), clonal_frequency=purity) 
    
    
    ##TIMER: VCF ##
    #solution from chatgpt to fix error with undefined superclass
    setOldClass("ExpData")
    
    vcf.snv <- VariantAnnotation::readVcf(vcfFileIn.snv, genome="hg38") 
    
    vcfFileIn.sv <- VariantAnnotation::readVcf(vcfFileIn.sv, genome="hg38")
    
    vcf.snv@info@listData$t_alt_count <- as.numeric(vcf.snv@info@listData$t_alt_count)
    vcf.snv@info@listData$t_ref_count <- as.numeric(vcf.snv@info@listData$t_ref_count)
    vcfFileIn.sv@info@listData$t_ref_count <- as.numeric(vcfFileIn.sv@info@listData$t_ref_count)
    vcfFileIn.sv@info@listData$t_alt_count <- as.numeric(vcfFileIn.sv@info@listData$t_alt_count)
    
    
    ssm.vcf <- VariantAnnotation::rbind(vcf.snv, vcfFileIn.sv)
    
    #### TIMER: run analysis #####
    cat('running core analysis... ')
    molecular_time_list <- mutationTime(ssm.vcf, cn, clusters=clusters, purity=purity, sex='female', isWgd=classWgd(cn), n.boot=10, xmin=3, rho=0)
    cat('done\n')
    #### TIMER: output ####
    
    vcf_df <- cbind.data.frame(rownames(info(ssm.vcf)), molecular_time_list$snv_timing, 'alt_count' =getAltCount(ssm.vcf), 'tumor_depth'=getTumorDepth(ssm.vcf))
    names(vcf_df)[1] <- 'ssm_id'
    type_list <- sapply(vcf_df, class)
    list_cols <- names(type_list[type_list == 'list'])
    vcf_df <- vcf_df[ , -which(names(vcf_df) %in% list_cols)]
    vcf_df$sample_id <- tumour_id
    
    write.table(
      vcf_df,
      file = paste0(working_dir,'/',tumour_id,'/mutationtimer/',tumour_id,'.sm.timed.txt'),
      row.names = FALSE, quote = FALSE, sep = "\t"
    )
    
    names(cnvs)[c(1:3)] <- c('chrom','start','end')
    cnv_df <- cbind.data.frame(cnvs[,c(1:3)], molecular_time_list$cn_timing)
    type_list <- sapply(cnv_df, class)
    list_cols <- names(type_list[type_list == 'list'])
    cnv_df <- cnv_df[ , -which(names(cnv_df) %in% list_cols)]
    
    write.table(
      cnv_df,
      file = paste0(working_dir,'/',tumour_id,'/mutationtimer/',tumour_id,'.cnv.timed.txt'),
      row.names = FALSE, quote = FALSE, sep = "\t"
    )
    
  } else {
    cat('some files missing for ',tumour_id,'\n')
    
  }
}

#### JOIN: make evolution file ####



first.clonality = T
for(tumour_id in sample_list$V2){
  
  pyclone.path <- paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/pyclone/',tumour_id,'.all.pyclone.tsv')
  depths.path <- paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/pyclone/',tumour_id,'.all.depths.txt')
  
  sigs.path <- paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/sigprofiler/',tumour_id,'.all.most_probable_signatures.txt')
  
  ## run with 0.1_run_mobster.R
  mobster.path <- paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/mobster/',tumour_id,'.all.mobster.txt')
  ## run with 0.2_clone_and_time.R
  timer.path <- paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/mutationtimer/',tumour_id,'.sm.timed.txt')
  
  if(all(file.exists(pyclone.path) &
         file.exists(depths.path) &
         file.exists(sigs.path) &
         file.exists(mobster.path) &
         file.exists(timer.path))
     
  ){
    
    dir_path = paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/evolution/')
    
    if (!dir.exists(dir_path)) {
      dir.create(dir_path)
    } 
    
    if(file.exists(paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/',tumour_id,'.evolution.txt'))){
      file.remove(paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/',tumour_id,'.evolution.txt'))
    }
    
    cat('processing ',tumour_id,'\n')
    
    pyclone <- fread(pyclone.path) %>% dplyr::select(mutation_id, cluster_id, cellular_prevalence, cluster_assignment_prob)
    mobster <- fread(mobster.path) %>% full_join(pyclone, by=c( 'mutation_id'='mutation_id'))
    names(mobster) <- c('sample_id','mutation_id','mobster.cluster','mobster.prob','pyclone.cluster','pyclone.prevalence','pyclone.prob')
    

    
    timer <- fread(timer.path) %>% full_join(mobster, by=c('sample_id'='sample_id','ssm_id'='mutation_id'))
    timer <- separate(timer, ssm_id, into = c('chr','pos','alt'), sep = ':', remove = F)
    timer$pos <- as.numeric(timer$pos)
    
    sigs <- fread(sigs.path) %>% full_join(timer, by=c('V1'='sample_id','V2'='chr', 'V3'='pos'))
    names(sigs)[1:6] <- c('tumour_id','chr','pos','context','sig','sig.prob')
    
    ##TODO: dnds_annotation does not contain any duplicate samples
    evolution.annotated <- sigs #left_join(sigs, dnds_annotation, by=c('tumour_id'='sampleID','chr'='chr','pos'='pos','alt'='mut'))
    
    evolution.annotated$chr <- factor(evolution.annotated$chr, levels = paste0('chr',c(c(1:22), c('X','Y'))))
    evolution.annotated <- evolution.annotated %>% arrange(chr, pos)
    
    write.table(
      evolution.annotated,
      file = paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/evolution/',tumour_id,'.evolution.raw.txt'),
      row.names = FALSE, quote = FALSE, sep = "\t"
    )
    
    chr.py <- evolution.annotated %>% group_by(chr, pyclone.cluster) %>% tally()
    
    py.clusters <- evolution.annotated %>% dplyr::select(pyclone.cluster, pyclone.prevalence)  %>% 
      dplyr::filter(!is.na(pyclone.cluster)) %>% unique() %>% arrange(desc(pyclone.prevalence))
    
    
    py.clusters$pyclone.rank <- seq(0,(nrow(py.clusters)-1))
    evolution.annotated <- left_join(evolution.annotated, py.clusters)
    evolution.annotated$fz <- evolution.annotated$alt_count /evolution.annotated$tumor_depth 
    
    evolution.annotated$sig[evolution.annotated$alt %in% sv.codes] <- 'sv'
    
    evolution.annotated$call <- NA
    if('Tail' %in% unique(evolution.annotated$mobster.cluster)){
      evolution.annotated$call[evolution.annotated$CLS == 'subclonal' & 
                                 evolution.annotated$mobster.cluster == 'Tail' &
                                 # evolution.annotated$pyclone.rank == max(py.clusters$pyclone.rank)
                                 evolution.annotated$pyclone.rank != 0 ] <- 'Tail'
      
      
    }
    
    evolution.annotated$call[evolution.annotated$CLS %in% clonal.times & 
                               evolution.annotated$mobster.cluster != 'Tail' &
                               evolution.annotated$pyclone.rank == 0 ] <- 'Clonal'
    
    
    evol.tally <- evolution.annotated %>% 
      dplyr::filter(chr %in% paste0('chr',c(1:22))) %>%
      group_by(tumour_id, CLS, mobster.cluster, pyclone.rank, call) %>%  #, pyclone.prevalence
      tally()
    
    evol.tally$percent <- evol.tally$n / sum(evol.tally$n)
    evol.tally$tumour_id[is.na(evol.tally$tumour_id)] <- tumour_id
    
    evol.filt <- evol.tally %>% dplyr::filter(percent > 0.02)
    
    significant.subclones <- evol.filt %>% 
      dplyr::filter(is.na(call) & 
               CLS == 'subclonal' & 
               mobster.cluster != 'Tail' &
               pyclone.rank != 0 )
    
    
    
    
    if(nrow(significant.subclones) > 0){
      if('Tail' %ni% unique(evolution.annotated$mobster.cluster)){
        ## makes the assumption that a tail is more likely than a single subclone when there is no tail
        
        evolution.annotated$call[is.na(evolution.annotated$call) & 
                                   evolution.annotated$CLS == 'subclonal' &
                                   evolution.annotated$mobster.cluster != 'Tail' &
                                   
                                   evolution.annotated$pyclone.rank == max(significant.subclones$pyclone.rank) ] <- 'Tail'
        
      } 
      
      if('Tail' %in% evol.filt$call){
        tail.clone <- evol.filt %>% dplyr::filter(call == 'Tail') %>% arrange(desc(percent) )
        
        tail.clone.id = tail.clone$pyclone.rank[1]
      } else {
        tail.clone.id = 'None'
      }
      
      for(this.subclone in unique(significant.subclones$pyclone.rank)){
        #ensures that disagreemens in the tail fit between mobster and pyclone
        #don't become a subclone
        if(this.subclone!=tail.clone.id){
          this.subclone.id = paste0('sub.',this.subclone)
          
          evolution.annotated$call[is.na(evolution.annotated$call) & 
                                     evolution.annotated$CLS == 'subclonal' &
                                     evolution.annotated$mobster.cluster != 'Tail' &
                                     
                                     evolution.annotated$pyclone.rank == this.subclone ] <- this.subclone.id
          
        }
      }
    }
    
    evol.tally.chr <- 
      evolution.annotated %>% 
      dplyr::filter(chr %in% paste0('chr',c(1:22))) %>%
      group_by(call, chr) %>%  
      summarize(count = n()) %>%
      group_by(call) %>%
      mutate(chr.fraction = count / sum(count))
    
    #evolution.annotated <- evolution.annotated %>% filter(!is.na(tumour_id))
    ggsave(paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/evolution/',tumour_id,'.evolution.raw.svg'),plot =
             
             print(
               cowplot::plot_grid(
                 
                 ggplot(evolution.annotated, aes(x=fz, fill=CLS)) + 
                   scale_fill_manual(values=c('lightblue','cornflowerblue','darkblue','forestgreen'), na.value = 'grey80') +
                   geom_histogram() + theme_bw() + theme(panel.grid = element_blank()) + labs(fill='TimeR') + xlim(0,1)
                 ,
                 ggplot(evolution.annotated, aes(x=fz, fill=mobster.cluster)) + geom_histogram()+ theme_bw()+ theme(panel.grid = element_blank()) + labs(fill='mobster')+ 
                   scale_fill_manual(values = pnw_palette(name="Bay", n=length(unique(evolution.annotated$mobster.cluster)) ), na.value = 'grey80')+ xlim(0,1)
                 ,
                 ggplot(evolution.annotated, aes(x=fz, fill=as.factor(pyclone.rank))) + geom_histogram()+ theme_bw()+ theme(panel.grid = element_blank()) + labs(fill='pyclone')+ 
                   scale_fill_manual(values = pnw_palette(name="Lake", n=length(unique(evolution.annotated$pyclone.rank)) ), na.value = 'grey80')+ xlim(0,1)
                 ,
                 ggplot(evolution.annotated, aes(x=fz, fill=call)) + geom_histogram()+ theme_bw()+ theme(panel.grid = element_blank()) + labs(fill='concensus')+ 
                   scale_fill_manual(values = pnw_palette(name="Sunset", n=length(unique(evolution.annotated$call)) ), na.value = 'grey80')+ xlim(0,1)
                 ,
                 ggplot(evol.tally.chr, aes(x=chr, y=call, fill=chr.fraction)) + geom_tile()+ theme_bw()+ theme(panel.grid = element_blank()) +
                   scale_fill_gradient2(  high = "cornflowerblue",
                                          low = "forestgreen", 
                                          mid = "white",
                                          na.value = "grey80") + labs(x='',y='',fill='')
                 ,
                 ncol = 1,align = 'hv', axis = "rlbt")
             )
           , device = "svg",width = 10, height = 10 )
    
    
    
    
    evolution.annotated$mut_length = nchar(evolution.annotated$alt)
    
    ggsave(paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/evolution/',tumour_id,'.evolution.sig.svg'),plot =
             
             print(
               
               ggplot(evolution.annotated %>% dplyr::filter(!is.na(call) & !is.na(sig) & mut_length == 1), aes(x=fz, fill=sig)) + 
                 geom_histogram()+ theme_bw(base_size=18)+ 
                 facet_grid(call~.)+
                 theme(panel.grid = element_blank()) + labs(fill='signature', x='Variant Allele Frequency', y='Mutation Count')+ 
                 scale_fill_manual(values = pnw_palette(name="Sunset", n=length(unique(evolution.annotated$sig[evolution.annotated$mut_length == 1])) ), na.value = 'grey80')+ xlim(0,1)
             )
           , device = "svg",width = 10, height = 6 )
    
    write.table(
      evolution.annotated,
      file = paste0('/Volumes/pcsi/users/fbeaudry/more.data/',tumour_id,'/evolution/',tumour_id,'.evolution.txt'),
      row.names = FALSE, quote = FALSE, sep = "\t"
    )
    
  } else {
    cat(tumour_id, ' not found\n')
    
  }
  

}
  
```  


