library(dplyr)
library(ACAT)
chr = 1
methods = c("lassosum", "P0.001", "P0.05", "PRScs", "SDPR", "FUSION")
res.dir = "~/projects/OTTERS/Data/ReRDA/Results"
anno.dir = "~/projects/OTTERS/Data/RDA/GTEx_V8/Anno"


ACAT_withNA = function(p_vec){
  p_vec_noNA = p_vec[is.na(p_vec) == F]
  ACAT(p_vec_noNA)
}

get_TWAS_p = function(chr){
  
  chr.anno.dir = file.path(anno.dir, paste0("GTEx_CHR", chr, "_GeneAnno.txt"))
  anno = read.table(chr.anno.dir, header = T)
  
  for (method in methods){
    method.res.dir = file.path(res.dir, 
                               paste0("CHR", chr),
                               paste0(method, ".txt"))
    res = read.table(method.res.dir, header = T, sep = '\t', row.names = NULL)  %>% 
      select(TargetID, FUSION_Z, FUSION_PVAL)
  
    #genomic control
    Lambda = median(res$FUSION_Z^2, na.rm = T)/qchisq(0.5, df = 1)
    res = res %>% mutate(p = pchisq(FUSION_Z^2/Lambda, df = 1, lower.tail = F)) %>% 
      select(TargetID, p)
    colnames(res) = c("TargetID", method)
    
    anno = merge(anno, res, by = "TargetID", all.x = T)
    
  }
  
  p_mat = anno[, c("lassosum", "P0.001", "P0.05", "PRScs", "SDPR")]
  p_acat = apply(p_mat, 1, ACAT_withNA)
  anno$ACAT = p_acat
  
  return(anno)
  
  
}

full.twas = file.path(res.dir, 'TWAS.txt')
chr1_twas = get_TWAS_p(chr = 1)

write.table(chr1_twas, full.twas, append = F,
            col.names = T, row.names = F, sep = "\t",
            quote = F)

for (chr in 2:22){
  
  chr_twas = get_TWAS_p(chr)
  
  write.table(chr_twas, full.twas, append = T,
              col.names = F, row.names = F, sep = "\t",
              quote = F)
  
}


