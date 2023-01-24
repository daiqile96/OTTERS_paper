#!/usr/bin/env Rscript
Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)

### Load in arguments
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  he2 = as.numeric(args[[1]])
  cp = as.numeric(args[[2]])
  N = as.numeric(args[[3]])
  hp2 = as.numeric(args[[4]])
  n_sim = as.numeric(args[[5]])
  n_set = as.numeric(args[[6]])
  n_thread = as.numeric(args[[7]])
  gene = as.character(args[[8]])
}

#----------------------SET UP-----------------------
.libPaths("/mnt/YangFSS/data/qdai/.local/R/4.0.4/lib64/R/library")
library(dplyr)
library(foreach)
library(bigstatsr)
library(MASS) 
library(doParallel)
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

get_Zscore <- function(j,causal_prop,heritability, N, hp2, seed = 2018){
  
  library(MASS)
  
  sub_res_path = paste0(out_path, "/std_cp_", causal_prop,"_He2_", heritability)
  true_beta_dir = paste0(sub_res_path,"/", gene, "_true_beta_", i, ".txt")
  true_beta_dat = read.table(true_beta_dir, header = T)
  true_beta = as.vector(true_beta_dat[, "BETA"])
  
  n_gwas = N * 1000
  mu = test_corr %*% true_beta * sqrt(n_gwas * hp2)
  
  Z = MASS::mvrnorm(n = 1, mu = mu, Sigma = test_corr)
  
  return(Z)
  
}

## Use parallelism by foreach in R
print(gene)

cl <- parallel::makeCluster(n_thread)
clusterEvalQ(cl, .libPaths("/mnt/YangFSS/data/qdai/.local/R/4.0.4/lib64/R/library"))
doParallel::registerDoParallel(cl)

foreach(k = 1:n_set, .combine = 'c') %dopar% {
  
  for (i in 1:n_sim){
    
    sub_res_path = paste0(out_path, "/std_cp_", cp,"_He2_", he2)
    
    Z = get_Zscore(i, cp, he2, N, hp2)
    
    snp_out$Z = as.vector(Z)
    
    # write Zscore
    write.table(snp_out,
                paste0(sub_res_path, "/", gene, "_Hp", hp2, "_Set", k, "_N", N, 
                       "K_GWAS_Zscore_",i, ".txt"),
                row.names=F, col.names=T, quote = F,sep = "\t")
  }
  
  print(paste("Finish simulate Zscore for sim", i))
  
}

parallel::stopCluster(cl)