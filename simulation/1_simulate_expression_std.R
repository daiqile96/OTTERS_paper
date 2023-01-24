#!/usr/bin/env Rscript

Sys.setlocale("LC_ALL", "C")
.libPaths("/mnt/YangFSS/data/qdai/.local/R/4.0.4/lib64/R/library")
library(parallel)
options(stringsAsFactors=F)
library(foreach)
library(bigstatsr)
library(doParallel)


### load in arguments
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  heritability = as.numeric(args[[1]])
  causal_prop = as.numeric(args[[2]])
  gene = as.character(args[[3]])
  n_thread = as.numeric(args[[4]])
  n_sim = as.numeric(args[[5]])
}

### settings
set.seed(8888)
data_path = "/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu"
out_path = file.path(data_path, "Training")
sub_out_path = paste0(out_path, "/std_cp_", causal_prop,"_He2_", heritability)
create_out_path = function(path){
  if (dir.exists(path) == FALSE){
    dir.create(path)
    return(paste0("Create ", path))
  } else {
    return(paste0(path, " exists"))
  }
}                 
lapply(c(out_path, sub_out_path), create_out_path)                  

# load the genotype data for gene 
geno_path = file.path(data_path, "binary", paste0("qc_", gene, "_geno.txt"))
geno_data = read.table(geno_path, header = T, check.names = F)

# extract snp information 
snp_data = geno_data[, 1:6] 
# extract genotype data
geno_data = geno_data[, 7:ncol(geno_data)] 
sample_id = colnames(geno_data)

# divided all samples into test (70%) and train data set (30%)
train_samples = as.vector(unlist(read.table(file.path(data_path, "train_id.txt"))))
test_samples = as.vector(unlist(read.table(file.path(data_path, "test_id.txt"))))
train_geno = geno_data[, train_samples]
test_geno = geno_data[, test_samples]
n_snp = nrow(geno_data)
n_sample = ncol(geno_data)
print("Complete loading data")


### Simulate Expression Data
print("Begin simulate expression data")
# select causal snps
n_causal = ceiling(n_snp * causal_prop)
print(paste0("n_causal = ", n_causal))
causal_snps <- sample(n_snp, n_causal)
# selected_snpinfo <- snp_info[selected_snps,]
std_geno_data = apply(geno_data, 1, scale, center = TRUE, scale = TRUE)
std_train_geno = apply(train_geno, 1, scale, center = TRUE, scale = TRUE)
causal_geno = as.matrix(std_geno_data[, causal_snps])

## Use parallelism by foreach in R
print("Begin parallel simulation")
cl <- parallel::makeCluster(n_thread)
clusterEvalQ(cl, .libPaths("/mnt/YangFSS/data/qdai/.local/R/4.0.4/lib64/R/library"))
doParallel::registerDoParallel(cl)
start_time <- Sys.time()

# simulate expression data & calculate summary statistics
foreach(i = 1:n_sim, .combine = 'c') %dopar% {
  
  library(dplyr)
  # generate a place holder 
  beta_matrix = as.matrix(rnorm(n_causal, mean = 0, sd = 1))
  
  # save true beta
  true_beta = rep(0, n_snp)
  true_beta[causal_snps] = as.vector(beta_matrix)
  
  # simulate expression data
  expr_initial <- causal_geno %*% beta_matrix
  # Calculate gamma 
  # VAR(gamma * X * Beta) = heritability
  gamma <- sqrt((heritability)/var(expr_initial[,1]))
  
  # the ture eQTL weights = gamma * true_beta
  snp_data$BETA = gamma * true_beta
  # write true eQTL weights
  write.table(snp_data, paste0(sub_out_path,"/", gene, "_true_beta_", i, ".txt"), 
              quote = F, row.names= F, col.names= T, sep = "\t")
  
  # Expr_causal_eQTL = gamma * (X*B)
  Expr1 <- gamma*expr_initial 
  # Error term = N(0, 1-heritability), so that VAR[gamma * X * Beta + error term] = 1
  stdev <- sqrt(1 - heritability)
  error_term <- rnorm(n_sample, mean = 0, sd = stdev) 
  # Overall Sim Expr = E1 + Err
  Expr = Expr1 + error_term
  Expr = data.frame(ID = sample_id, Expr = Expr)
  rownames(Expr) = sample_id
  
  # get expression data for training and test samples
  train_Expr = Expr[train_samples,]
  test_Expr = Expr[test_samples,]
  
  cal_sst <- function(exp,geno,i){
    
    snp_vec <- unlist(geno[,i])
    exp_vec <- unlist(exp)
    
    res = summary(lm(exp_vec~snp_vec))
    # beta = res$coefficients[2, 1]
    # p = res$coefficients[2, 4]
    z = res$coefficients[2, 3]
    
    return(z)
  }
  
  # standardize expression data for training samples
  std_train_exp = scale(unlist(train_Expr[,"Expr"]), center = T, scale = T)
  # calculate eQTL summary statistics using training samples
  eQTL_Z = sapply(1:n_snp, function(i) cal_sst(std_train_exp, std_train_geno, i))
  
  # save eQTL summary statistics (Z score) into a data frame
  colnames(snp_data) = c("CHROM", "SNP", "POS", "A1", "A2", "MAF", "True_BETA")
  snp_data = snp_data %>% 
    mutate(Z = eQTL_Z,
           TargetID = gene,
           N = nrow(std_train_geno)) %>% 
    select(CHROM, POS, A1, A2, Z, TargetID, N)
  
  # write simulated expression data for test samples
  write.table(test_Expr, paste0(sub_out_path, "/", gene, "_test_Expr_", i, ".txt"),
              quote = F, row.names = F, col.names = T, sep = "\t")

  # save simulated expression data for test samples
  write.table(train_Expr, paste0(sub_out_path, "/", gene, "_train_Expr_", i, ".txt"),
              quote = F, row.names= F, col.names= T, sep = "\t")

  # write eQTL summary statistics
  write.table(snp_data, paste0(sub_out_path,"/", gene, "_sst_", i, ".txt"),
              quote = F, row.names= F, col.names= T, sep = "\t")
  
}

end_time <- Sys.time()
# print(paste("Computation time: ", end_time - start_time))

parallel::stopCluster(cl)