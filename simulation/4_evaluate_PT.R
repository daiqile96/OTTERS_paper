#!/usr/bin/env Rscript
.libPaths("/mnt/YangFSS/data/qdai/.local/R/4.0.4/lib64/R/library/")
Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)
library(parallel)
library(dplyr)
library(ACAT)

### Load in arguments
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  he2 = as.numeric(args[[1]])
  cp = as.numeric(args[[2]])
  N = as.numeric(args[[3]])
}

##################### Get group information ########################
Simu.dir = "~/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu"

# Read in gene annotation file (in hg19)
GeneLength = read.table(file.path(Simu.dir, "0_GeneLength.txt"), header = T)
GeneLength = GeneLength %>% 
  filter(CHROM %in% c(1:22)) %>% 
  mutate(Length = GeneEnd - GeneStart) %>% 
  arrange(Length) 

# select by quantile
length.qtl = quantile(GeneLength$Length, probs = seq(0, 1, 0.2))

GeneLength = GeneLength %>% 
  mutate(length.group = case_when(
    Length < 9368.4 ~ 1,
    Length >= 9368.4 & Length < 22368.8 ~ 2,
    Length >= 22368.8 & Length < 45435.6 ~ 3,
    Length >= 45435.6 & Length < 102263.6 ~ 4,
    Length >= 102263.6 & Length < 5379014.0 ~ 5
  ))

# merge with the selected genes
gene.anno = read.table(file.path(Simu.dir, "0_ReSimu_GeneAnno.txt"), header = F)
colnames(gene.anno) = c("CHROM", "GeneStart", "GeneEnd", "TargetID", "GeneName", "Length")
gene.anno = merge(gene.anno, GeneLength[, c("TargetID", "length.group")], by = "TargetID")

##################### Extract power ########################

get_Power = function(gene){
  
  get_power_method_set = function(method, set){
    
    # read in estimated GReX
    power.dir = file.path(TWAS.dir, paste0(gene, "_Set", set), paste0(method, '.txt'))
    
    tmp.p = NA
    
    if (file.exists(power.dir)){
      
      power.tab = read.table(power.dir, header = T, check.names = F, stringsAsFactors = F, fill = T) 
      
      if (nrow(power.tab) > 0 ){
        
        if (is.na(power.tab[, 'FUSION_Z']) == F){
          
          tmp.p = 2*pnorm(q=abs(power.tab[, 'FUSION_Z']), lower.tail=FALSE)
          
        }
        
      } 
      
    }
    
    return(tmp.p)
    
  }
  
  get_power_set = function(set){
    
    tmp.res = unlist(sapply(methods, get_power_method_set, set = set))
    
    prscs.res = as.numeric(tmp.res['PRScs'])
    pt0.01.res = as.numeric(tmp.res['P0.001'])
    
    if (is.na(pt0.01.res) & is.na(prscs.res) == F){
      
      idx = case_when(
        prscs.res < 2.5e-6 ~ 4,
        prscs.res >= 2.5e-6 ~ 5,
      )
      
      
    } else {
      
      idx = case_when(
        prscs.res < 2.5e-6 & pt0.01.res < 2.5e-6 ~ 1,
        prscs.res >= 2.5e-6 & pt0.01.res >= 2.5e-6 ~ 0,
        prscs.res < 2.5e-6 & pt0.01.res >= 2.5e-6 ~ 2,
        prscs.res >= 2.5e-6 & pt0.01.res < 2.5e-6 ~ 3,
      )
      
      return(idx)
      
    }
    
  }
  
  res = sapply(1:10, get_power_set)
  
  get_smallest_p_method = function(p_vec){
    
    return(names(which.min(p_vec)))
    
  }
  
  # power.res = rowSums(res < 2.5e-6, na.rm = T)
  print(gene)
  
  return(res)
  
  
}

res_by_group = function(group){
  
  genes = gene.anno %>% filter(length.group == group) %>% select(GeneName) %>% unlist() %>% unname()
  res = sapply(genes, get_Power) 
  
  res.power = rowSums(res)/1000
  names(res.power) = c(methods, "ACAT")
  res.vec = c(cp, he2, N, group, unlist(res.power))
  
  write.table(t(res.vec), 
              file.path(Simu.dir, "Results", "0_ReSimu_Power_by_group.txt"),
              row.names = F,
              col.names = F,
              quote = F,
              sep = '\t',
              append = T)
}

sub.dir = paste0("std_cp_", cp, "_He2_", he2)
TWAS.dir = file.path(Simu.dir, "Results", sub.dir, "TWAS", paste0("N", N))
methods = c("SDPR", "PRScs", "lassosum", "P0.05", "P0.001")
cmp.res = sapply(gene.anno$GeneName, get_Power)

table(cmp.res)

