#!/usr/bin/env Rscript

Sys.setlocale("LC_ALL", "C")
options(stringsAsFactors=F)

### load in arguments
args=(commandArgs(TRUE))
print(args)
if(length(args)==0) {
  stop("Error: No arguments supplied!")
} else {
  chr = as.numeric(args[[1]])
  batch = as.numeric(args[[2]])
}

library(dplyr)
project.dir = "~/projects/OTTERS/Data"
bim.dir = file.path(project.dir, "ReRDA/ReFUSION", 
                    paste0("CHR", chr))
anno.dir = file.path(project.dir, "RDA/GTEx_V8/Anno/", 
                     paste0("Batch", batch, "_GTEx_CHR", chr , "_GeneAnno.txt"))
anno = read.table(anno.dir, header = T)
out.dir = file.path(project.dir, "ReRDA/Results", 
                    paste0("CHR", chr), paste0("Batch", batch))

# Create output file
out_cols = c("CHROM", "POS", "A1", "A2", "TargetID", "ES")
MyDf = data.frame(matrix(nrow = 0, ncol = length(out_cols)))
colnames(MyDf) = out_cols

out_file = file.path(out.dir, paste0("FUSION", ".txt"))
write.table(MyDf, file = out_file,
            row.names = F,
            col.names = T,
            append = F,
            quote = F,
            sep = "\t")

for (ID in anno$TargetID){
  
  fusion.out = file.path(bim.dir, ID, paste0(ID, ".wgt.RDat"))
  
  if (file.exists(fusion.out)){

    #cat( unlist(wgtlist[w,]) , '\n' )
    # Load weights
    load(fusion.out)
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
      cat( "WARNING : " , unlist(wgtlist[w,]) , " did not have a predictive model ... skipping entirely\n" )
      FAIL.ctr = FAIL.ctr + 1
      next
    }
    
    if ( names(mod.best) == "lasso" || names(mod.best) == "enet" ) {
      keep = wgt.matrix[,mod.best] != 0
    } else if ( names(mod.best) == "top1" ) {
      keep = which.max(wgt.matrix[,mod.best]^2)
    } else { 
      keep = 1:nrow(wgt.matrix)
    }
    
    if (sum(keep) == 0){
      
      cat("Target", ID, "No Fusion Results.\n", sep = " ")
      
    } else {
      
      # Create a data frame including SNPs and estimated eQTL weights
      wgt.df = data.frame(SNP = row.names(wgt.matrix)[keep], 
                          ES = wgt.matrix[keep, mod.best],
                          TargetID = ID)
      
      # Format snps file 
      colnames(snps) = c("CHROM", "SNP", "BP", "POS", "A1", "A2")
      
      # Merge snps with wgt.df to get estimated eQTL weights for each SNP
      df_out = merge(snps, wgt.df, by = "SNP")
      
      df_out = df_out %>% arrange(as.numeric(POS))
      
      # Save the weights
      write.table(df_out[, out_cols], 
                  file = out_file, 
                  row.names = F,
                  col.names = F,
                  append = T,
                  quote = F,
                  sep = "\t")
      
      
      print(paste("Done Format Fusion", ID))
      
      
    }
    
  } else {
    
    cat("Target", ID, "No Fusion Results.\n", sep = " ")
    
  }
  
}

