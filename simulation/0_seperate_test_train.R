# This script seperates the 1894 samples into 568 training samples and 1326 test samples
# seed: 8888
# output: train_id.txt and test_id.txt including the sample ID of training samples and test samples

# read in genotype data for gene ABCA7
data_path = "/mnt/YangFSS/data2/qdai/BS_TWAS/Data/Simulation"
gene <- "ABCA7"

# load the genotype data for gene ABCA7 
geno_path = file.path(data_path, paste0(gene, "_genotype.txt"))
geno_data = read.table(geno_path, header = T, check.names = F)

# extract genotype data
geno_data = geno_data[, 7:ncol(geno_data)] 
sample_id = colnames(geno_data)

# divided all samples into test (70%) and train data set (30%)
set.seed(8888)
train_size = floor(length(sample_id) * 0.3)
train_samples = sample(sample_id, size = train_size, replace = F)
test_samples = sample_id[!sample_id %in% train_samples]

# use new out path
out_path = "/home/qdai/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu"
# save sample id of test samples and training samples
write.table(train_samples, file.path(out_path, "train_id.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
write.table(test_samples, file.path(out_path, "test_id.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
