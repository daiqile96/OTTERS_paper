#!/usr/bin/env Rscript

Sys.setlocale("LC_ALL", "C")
.libPaths("/home/qdai8/projects/Rlib")
options(stringsAsFactors=F)

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

# Fusion Results
res.dir = "~/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu/Results/"
FUSION.dir = file.path(res.dir, "FUSION")
FUSION.res = file.path(FUSION.dir, paste0(gene, "_", cp, "_", he2, ".wgt.RDat"))

print(FUSION.res)

# Create output file
out_cols = c("CHROM", "POS", "A1", "A2", "TargetID", "ES")
MyDf = data.frame(matrix(nrow = 0, ncol = length(out_cols))) 
colnames(MyDf) = out_cols


## For each wgt file:
if (file.exists(FUSION.res)){
  
  # Create output file
  sub.dir = file.path(res.dir, paste0('std_cp_', cp, '_He2_', he2), gene)
  out_file = file.path(sub.dir, paste0("FUSION.txt"))
  
  write.table(MyDf, file = out_file,
              row.names = F,
              col.names = T,
              append = F,
              quote = F,
              sep = "\t")
  
  #cat( unlist(wgtlist[w,]) , '\n' )
  # Load weights
  load(FUSION.res)
  # Remove NAs (these should not be here)
  wgt.matrix[is.na(wgt.matrix)] = 0
  
  # which rows have rsq
  row.rsq = grep( "rsq" , rownames(cv.performance) )
  # which rows have p-values
  row.pval = grep( "pval" , rownames(cv.performance) )	
  
  # Identify the best model
  # get the most significant model
  mod.best = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))
  
  print(cv.performance[row.pval,])
  
  if ( length(mod.best) == 0 ) {
    cat( "WARNING : " , gene, " did not have a predictive model ... skipping entirely\n" )
    FAIL.ctr = FAIL.ctr + 1
    next
  }
  
  # Create a data frame including SNPs and estimated eQTL weights
  wgt.df = data.frame(SNP = row.names(wgt.matrix), ES = wgt.matrix[, mod.best],
                      TargetID = gene)
  
  # Format snps file 
  colnames(snps) = c("CHROM", "SNP", "BP", "POS", "A1", "A2")
  
  # Merge snps with wgt.df to get estimated eQTL weights for each SNP
  df_out = merge(snps, wgt.df, by = "SNP")
  
  # Save the weights
  write.table(df_out[, out_cols], 
              file = out_file, 
              row.names = F,
              col.names = F,
              append = T,
              quote = F,
              sep = "\t")
  
  print(paste("Done res", gene))
  
} else {
  
  print(paste(FUSION.res, "doesn't exist"))
  
}
