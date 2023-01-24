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
# gene.anno = read.table(file.path(Simu.dir, "0_ReSimu_GeneAnno.txt"), header = F)
gene.anno = read.table(file.path(Simu.dir, "0_ReSimu_GeneAnno_100_per_group.txt"), header = F)
colnames(gene.anno) = c("CHROM", "GeneStart", "GeneEnd", "TargetID", "GeneName", "Length")
gene.anno = merge(gene.anno, GeneLength[, c("TargetID", "length.group")], by = "TargetID")

##################### Extract power ########################
ACAT_with_NA = function(p_vec){
  
  p_vec = p_vec[is.na(p_vec) == F]
  
  return(ACAT(p_vec))
  
}

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
    
    # no FUSION 
    ACAT.res = ACAT_with_NA(tmp.res[1:5])
      
    tmp.res = c(tmp.res, ACAT.res)
    
    return(tmp.res)
  }
  
  res = sapply(1:10, get_power_set)
  res = apply(res, 2, as.numeric)
  
  power.res = rowSums(res < 2.5e-6, na.rm = T)
  
  return(power.res)
  
  
}

##################### Extract power by group ########################
power_by_group = function(group){
  genes = gene.anno %>% filter(length.group == group) %>% select(GeneName) %>% unlist() %>% unname()
  res = sapply(genes, get_Power) 
  res.power = rowSums(res)/(length(genes) * 10)
  names(res.power) = c(methods, "ACAT")
  res.vec = c(cp, he2, N, group, unlist(res.power))
  
  write.table(t(res.vec), 
              file.path(Simu.dir, "Results", "2_ReSimu_Power_by_group_100_per_group.txt"),
              row.names = F,
              col.names = F,
              quote = F,
              sep = '\t',
              append = T)
}

sub.dir = paste0("std_cp_", cp, "_He2_", he2)
TWAS.dir = file.path(Simu.dir, "Results", sub.dir, "TWAS", paste0("N", N))
methods = c("SDPR", "PRScs", "lassosum", "P0.05", "P0.001", "FUSION")
sapply(1:5, power_by_group)

# for all groups
genes = gene.anno %>% select(GeneName) %>% unlist() %>% unname()
res = sapply(genes, get_Power) 
res.power = rowSums(res)/(length(genes) * 10)
names(res.power) = c(methods, "ACAT")
res.vec = c(cp, he2, N, unlist(res.power))
write.table(t(res.vec), 
            file.path(Simu.dir, "Results", "2_ReSimu_Power_100_per_group.txt"),
            row.names = F,
            col.names = F,
            quote = F,
            sep = '\t',
            append = T)
