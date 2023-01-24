#!/usr/bin/env Rscript

Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)
suppressMessages(library("dplyr"))

### load in arguments
args=(commandArgs(TRUE))
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  chr = as.numeric(args[[1]])
}

# # Save expression data with clean TargetID
# Exp = read.table('~/projects/OTTERS/Data/RDA/GTEx_V6/Whole_Blood_Analysis.v6p.normalized.expression.bed',
#                  comment.char = "!",
#                  check.names = F,
#                  header = T)
# colnames(Exp)[1:4] = c("CHROM", "GeneStart", "GeneEnd", "TargetID")
# cleanid <- function(gtex_id) sapply(strsplit(gtex_id, "\\."), "[", 1L)
# Exp$TargetID = cleanid(Exp$TargetID)
# write.table(Exp, file.path('/home/qdai8/projects/OTTERS/Data/RDA/GTEx_V6',
#                            'Whole_Blood_GTEx_Exp.txt'),
#             sep = '\t',
#             col.names = T,
#             row.names = F,
#             quote = F)

# settings
chr.str = paste0("CHR", chr)
Data.dir = "~/projects/OTTERS/Data"
bim.dir = file.path(Data.dir, "RDA/FUSION", chr.str)
anno.dir = file.path(Data.dir, "RDA/GTEx_V8/Anno", paste('GTEx', chr.str, 'GeneAnno.txt',
                                                     sep = '_'))
anno = read.table(anno.dir, header = T)

# read in true expression
exp.dir = file.path(Data.dir, "RDA/GTEx_V6", 'Whole_Blood_GTEx_Exp.txt')
Exp = read.table(exp.dir, sep = '\t', header = T, check.names = F)

# read in ID of training samples
Train.ID.dir = file.path(Data.dir, "RDA/GTEx_V6", 'TrainID.txt')
Train.ID = read.table(Train.ID.dir)[,1] %>% unlist()
rownames(Exp) = Exp$TargetID

# change expression data in fam file
change_fam_exp = function(gene){
  
  fam.dir = file.path(bim.dir, gene, paste0(gene, ".fam"))
  fam = data.frame(read.table(fam.dir, header = F))
  colnames(fam) = c("FID", "ID", "ID_F", "ID_M", "SEX", "Pheno")
  
  exp.vec = unlist(Exp[gene, Train.ID])
  exp = data.frame(Expr = exp.vec, ID = Train.ID)
  
  # generate new fam file (with expression data)
  new_fam = merge(fam, exp, by = "ID")
  
  if (nrow(new_fam) != nrow(fam)){
    print('Unequal samples')
  }
  
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

sapply(anno$TargetID, change_fam_exp)