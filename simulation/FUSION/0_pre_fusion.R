#!/usr/bin/env Rscript

Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)
suppressMessages(library("dplyr"))

### load in arguments
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  cp = as.numeric(args[[1]])
  he2 = as.numeric(args[[2]])
  gene = as.character(args[[3]])
}

bim.dir = "~/projects/OTTERS/Data/ReSimu"

# paths to expression data
sub.dir = file.path(bim.dir, 'Training', paste0('std_cp_', cp, '_He2_', he2))


# change expression data in fam file
change_fam_exp = function(gene){
  
  fam.dir = file.path(bim.dir, "binary", paste0("train_", gene, ".fam"))
  fam = data.frame(read.table(fam.dir, header = F))
  colnames(fam) = c("FID", "ID", "ID_F", "ID_M", "SEX", "Pheno")
  
  exp.dir = file.path(sub.dir, paste0(gene, '_train_Expr_1.txt'))
  exp = data.frame(read.table(exp.dir, header = T))
  
  # generate new fam file (with expression data)
  new_fam = merge(fam, exp, by = "ID")
  row.names(new_fam) = new_fam$ID
  new_fam = new_fam[fam$ID, ] %>% 
    select(FID, ID, ID_F, ID_M, SEX, Expr)
  
  write.table(new_fam,
              fam.dir,
              col.names = F,
              quote = F,
              sep = "\t",
              row.names = F)
  
}

change_fam_exp(gene)

print('Finish Generating Fam with Expression Data.')
