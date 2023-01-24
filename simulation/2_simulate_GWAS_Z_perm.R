#!/usr/bin/env Rscript
Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)

### Load in arguments
args=(commandArgs(TRUE))
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  he2 = as.numeric(args[[1]])
  cp = as.numeric(args[[2]])
  n_sim = as.numeric(args[[3]])
  n_set = as.numeric(args[[4]])
  n_thread = as.numeric(args[[5]])
  gene = as.character(args[[6]])
}

#----------------------SET UP-----------------------
.libPaths("/mnt/YangFSS/data/qdai/.local/R/4.0.4/lib64/R/library")
set.seed(2018)
library(dplyr, quietly = T, warn.conflicts = F)
library(foreach, quietly = T, warn.conflicts = F)
library(bigstatsr, quietly = T, warn.conflicts = F)
library(MASS, quietly = T, warn.conflicts = F) 
library(doParallel, quietly = T, warn.conflicts = F)
data_path = "/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu"

#----------------------Read In Genotype Data & Sample ID-----------------------

### settings
out_path = file.path(data_path, "Training")

# load the genotype data for gene ABCA7 
geno_path = file.path(data_path, "binary", paste0("qc_", gene, "_geno.txt"))
geno_data = read.table(geno_path, header = T, check.names = F)

# extract snp information 
snp_data = geno_data[, 1:6] 
# extract genotype data
geno_data = geno_data[, 7:ncol(geno_data)] 
sample_id = colnames(geno_data)

# divided all samples into test (70%) and train data set (30%)
test_samples = as.vector(unlist(read.table(file.path(data_path, "test_id.txt"))))
test_geno = geno_data[, test_samples]
n_snp = nrow(geno_data)

# Standardize genotype data 
std_test_geno = apply(test_geno, 1, scale, center = TRUE, scale = TRUE)
colnames(std_test_geno) = snp_data[, "SNP"]
rownames(std_test_geno) = test_samples
n_test = length(test_samples)
# std_test_geno (1326 2156)

# calculate covariance matrix in test samples
test_corr = cov(std_test_geno)

# extract SNP information
snp_out = snp_data[, c("Chrom", "SNPPos", "A1", "A2")]
colnames(snp_out) = c("CHROM", "POS", "A1", "A2")

#-----------------------Functions-----------------------------

cal_Zscore <- function(i,pheno){
  #----------------------------------------------------
  # read in and scale the true genotype data for one SNP
  #----------------------------------------------------
  std_geno = as.vector(unlist(std_test_geno[,i]))
  
  #----------------------------------------------------
  # calculate GWAS z-score for on SNP
  #----------------------------------------------------
  res = summary(lm(pheno~std_geno))
  Z = res$coefficients[2, 3]
  return(Z)
  
}

get_Zscore <- function(j){
  
  test_pheno = rnorm(n_test)
  std_test_pheno = scale(test_pheno, center = TRUE, scale = TRUE)
  Z = sapply(1:n_snp, function(i) cal_Zscore(i, std_test_pheno))
  
  return(Z)
  
}

## Use parallelism by foreach in R
print(gene)


for (i in 1:n_sim){
  
  sub_res_path = paste0(out_path, "/std_cp_", cp,"_He2_", he2)
  
  Z = get_Zscore(i)
  
  snp_out$Z = as.vector(Z)
  
  write.table(snp_out,
              paste0(sub_res_path, "/Null_", gene, "_Set", 1, 
                     "_GWAS_Zscore_",i, ".txt"),
              row.names=F, col.names=T, quote = F,sep = "\t")
  
}

