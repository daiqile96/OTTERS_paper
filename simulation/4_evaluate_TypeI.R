library(ACAT)
library(dplyr)

Simu.dir = "~/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu"
# gene.anno = read.table(file.path(Simu.dir, "0_ReSimu_GeneAnno.txt"), header = F)
gene.anno = read.table(file.path(Simu.dir, "0_ReSimu_GeneAnno_100_per_group.txt"), header = F)
colnames(gene.anno) = c("CHROM", "GeneStart", "GeneEnd", "TargetID", "GeneName", "Length")
genes = gene.anno %>% select(GeneName) %>% unlist() %>% unname()

# permutation
Perm.dir = file.path(Simu.dir, '/Results/std_cp_0.001_He2_0.1/TWAS/Perm')
methods = c("SDPR", "PRScs", "lassosum", "P0.05", "P0.001")

ACAT_with_NA = function(p_vec){
  
  p_vec = p_vec[is.na(p_vec) == F]
  
  return(ACAT(p_vec))
  
}

get_p = function(gene){
  
  get_p_method = function(method){
    
    # read in estimated GReX
    power.dir = file.path(Perm.dir, gene, paste0(method, '.txt'))
    
    tmp.p = rep(NA, 2000)
    
    if (file.exists(power.dir)){
      
      power.tab = read.table(power.dir, header = T, check.names = F, stringsAsFactors = F, fill = T) 
      
      if (nrow(power.tab) == 2000){
        
        if (sum(is.na(power.tab[, 'FUSION_Z'])) == 0){
          
          tmp.p = 2*pnorm(q=abs(power.tab[, 'FUSION_Z']), lower.tail=FALSE)
          
        }
        
      } 
      
    }
    
    return(tmp.p)
    
  }
  
  tmp.res = unlist(sapply(methods, get_p_method))
  ACAT.res = apply(tmp.res, 1, ACAT_with_NA)
  tmp.res = cbind(tmp.res, ACAT.res)
  
  return(tmp.res)
  
  
}

res = lapply(genes, get_p)
saveRDS(res, file = file.path(Perm.dir, 'Perm_res.rds'))

