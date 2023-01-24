#!/usr/bin/env Rscript

Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)

### load in arguments
args=(commandArgs(TRUE))
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  twas.res.dir = as.character(args[[1]])
  GReX.r2.dir = as.character(args[[2]])
  anno.dir = as.character(args[[3]])
  methods = strsplit(args[[4]], ",")[[1]]
  filter = as.character(args[[5]])
}

# methods = c("lassosum", "P0.001", "P0.05", "PRScs", "SDPR")
# twas.res.dir = "~/projects/OTTERS/Data/RDA/Results/TWAS"
# GReX.r2.dir = "~/projects/OTTERS/Data/RDA/Results/GReX_R2.txt"
# anno.dir = "~/projects/OTTERS/Data/RDA/GTEx_V8/Anno"


if (filter == "T"){
  out.twas.dir = file.path(twas.res.dir, 'TWAS_Filter.txt')
} else {
  out.twas.dir = file.path(twas.res.dir, 'TWAS.txt')
}


######################## SET UP ################################
library(dplyr)
library(ACAT)

######################## FUNCTIONS #############################
# Perform ACAT with NA in p-values
ACAT_withNA = function(p_vec){
  p_vec_noNA = p_vec[is.na(p_vec) == F]
  ACAT(p_vec_noNA)
}

# Get TWAS p-values for methods (filtered results with test R2 < test_R2_cutoff)
get_filtered_TWAS = function(chr, test_R2_cutoff, methods){
  
  # read in annotation file 
  chr.anno.dir = file.path(anno.dir, paste0("GTEx_CHR", chr, "_GeneAnno.txt"))
  anno = read.table(chr.anno.dir, header = T)
  
  for (method in methods){
    
    ################# read and format the TWAS results ########
    method.res.dir = file.path(twas.res.dir, 
                               paste0("CHR", chr),
                               paste0(method, ".txt"))
    res = read.table(method.res.dir, header = T, sep = '\t', row.names = NULL)  %>% 
      select(TargetID, FUSION_Z, FUSION_PVAL)
    
    ################# genomic control #########################
    Lambda = median(res$FUSION_Z^2, na.rm = T)/qchisq(0.5, df = 1)
    res = res %>% mutate(p = pchisq(FUSION_Z^2/Lambda, df = 1, lower.tail = F)) %>% 
      select(TargetID, p)
    colnames(res) = c("TargetID", method)
    
    ################# merge with annotation ###################
    anno = merge(anno, res, by = "TargetID", all.x = T)
    
  }
  
  ################# get filterd p-values ##################
  # get the TWAS p-value for all the methods
  p_mat = anno[, methods]
  rownames(p_mat) = anno$TargetID
  
  if (filter == "T"){
    # get the test R2 matrix
    sub.r_mat = r_mat[rownames(p_mat), ]
    # if test R2 >= test_R2_cutoff, set as 1, else NA
    idx = sub.r_mat < test_R2_cutoff
    sub.r_mat[idx] = NA
    sub.r_mat[!idx] = 1
    # get TWAS p-value for all the methods after filter results with test R2 < 0.01
    sub.p_mat = p_mat * sub.r_mat
    
    ################## ACAT ######################
    ACAT_methods = c('P0.001', 'P0.05', 'lassosum', 'SDPR', 'PRScs')
    acat.p_mat = sub.p_mat[, ACAT_methods]
    
    ################## replace the raw TWAS results with filtered results #######
    anno[, methods] = sub.p_mat
    
  } else {
    
    ################## ACAT ######################
    ACAT_methods = c('P0.001', 'P0.05', 'lassosum', 'SDPR', 'PRScs')
    acat.p_mat = p_mat[, ACAT_methods]
    
  }
  
  p_acat = apply(acat.p_mat, 1, ACAT_withNA)
  anno$ACAT = p_acat
  
  anno
  
}

######################## GET TWAS RESULTS #############################
# Extract test R2 matrix
GReX.R2 = read.table(GReX.r2.dir, 
                     header = T, sep = '\t')
r_mat = GReX.R2[, methods]
rownames(r_mat) = GReX.R2$TargetID

# Output Results
for (chr in 1:22){
  
  chr_twas = get_filtered_TWAS(chr = chr, test_R2_cutoff = 0.01, methods)
  
  if (chr == 1){
    write.table(chr_twas, out.twas.dir, append = F,
                col.names = T, row.names = F, sep = "\t",
                quote = F)
  } else{
    write.table(chr_twas, out.twas.dir, append = T,
                col.names = F, row.names = F, sep = "\t",
                quote = F)
  }

}

