#! /usr/bin/env Rscript

reformat_cox <- function(palliative.cohort, cox.os.cms.loc.palliative){
  
  cox.os.cms.loc.palliative.summary <-   summary(cox.os.cms.loc.palliative)
  
  full.forest.df <- cbind.data.frame(
    cox.os.cms.loc.palliative.summary$conf.int,
    'p'=cox.os.cms.loc.palliative.summary$coefficients[,5]
  ) %>% as.data.frame() %>% 
    dplyr::select(`exp(coef)`, `lower .95`, `upper .95`, p)
  
  names(full.forest.df) <- c('HZ','ci.low','ci.hi','p')
  full.forest.df$name <- rownames(full.forest.df)
  
  palliative.n <- palliative.cohort %>% 
    dplyr::select(rna_class, location, Overall_Stage_at_Diagnosis ) %>% 
    melt(measure.vars=c('rna_class', 'location', 'Overall_Stage_at_Diagnosis'))
  
  palliative.n$name <- paste0(palliative.n$variable, palliative.n$value )
  palliative.n <- palliative.n %>% group_by(name, variable, value)  %>% 
    filter(!is.na(value)) %>% tally()
  
  full.forest.df <- left_join(palliative.n, full.forest.df)
  full.forest.df$HZ[is.na(full.forest.df$HZ)] <- 1
  
  full.forest.df$p <- round(full.forest.df$p,3)
  
  full.forest.df <- full.forest.df %>%
    mutate(p = if_else(!is.na(p) & p < 0.05,
                       paste0(p, "*"),
                       as.character(p)))
  
  full.forest.df$p[is.na(full.forest.df$p)] <- '(reference)'
  
  full.forest.df$variable <- as.character(full.forest.df$variable)
  full.forest.df$variable[full.forest.df$variable == 'Overall_Stage_at_Diagnosis'] <- 'Stage'
  full.forest.df$variable[full.forest.df$variable == 'rna_class'] <- 'eCCS'
  full.forest.df$category <- factor(full.forest.df$variable, levels = c('eCCS','location','Stage'))
  
  full.forest.df$value[full.forest.df$value == 'Distal'] <- 'dCCA'
  full.forest.df$value[full.forest.df$value == 'Gallbladder'] <- 'GBC'
  full.forest.df$value[full.forest.df$value == 'Intrahepatic'] <- 'iCCA'
  full.forest.df$value[full.forest.df$value == 'Perihilar'] <- 'pCCA'
  
  full.forest.df$value.n <- paste0(full.forest.df$value,'\nn=',full.forest.df$n)
  
  full.forest.df$value <- factor(full.forest.df$value, levels = c('IV','III','I & II',
                                                                  'dCCA','GBC','pCCA','HCC-CCA','iCCA',
                                                                  'CCS-B','CCS-A'))
  
  full.forest.df$type <- factor(full.forest.df$value.n, levels = full.forest.df$value.n[order(full.forest.df$value)])
  return(full.forest.df)
}
