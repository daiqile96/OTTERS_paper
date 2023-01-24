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
}

Simu.dir = "~/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu"
sub.dir = paste0("std_cp_", cp, "_He2_", he2)
GReX.dir = file.path(Simu.dir, "Results", sub.dir, "GReX")
True.dir = file.path(Simu.dir, "Training",  sub.dir)

gene.anno = read.table(file.path(Simu.dir, "0_ReSimu_GeneAnno_100_per_group.txt"), header = F)
genes = unlist(gene.anno$V5)

methods = c("SDPR", "PRScs", "lassosum", "P0.05", "P0.001", "FUSION")

get_R2 = function(gene){
  
  trueGReX.dir = file.path(True.dir, paste0(gene, "_test_Expr_1.txt"))
  
  if (file.exists(trueGReX.dir)){
    trueGReX = read.table(trueGReX.dir, header = T)
    
    get_R2_one_method = function(method){
      
      # read in estimated GReX
      estGReX.dir = file.path(GReX.dir, gene, paste0(method, ".txt"))
      
      
      if (file.exists(estGReX.dir)){
        
        estGReX = data.frame(t(read.table(estGReX.dir, header = T, check.names = F, row.names = 1,
                                          stringsAsFactors = F))) 
        
        if (ncol(estGReX) > 1){
          # format estimated GReX
          colnames(estGReX) = c("IID", "EstExpr")
          estGReX$ID = unlist(row.names(estGReX))
          estGReX$EstExpr = as.numeric(estGReX$EstExpr)
          
          # merge with true expression data
          allGReX = merge(estGReX, trueGReX, by = "ID")
          
          # calculate R2
          res = cor.test(allGReX$EstExpr, allGReX$Expr)
          return(res$estimate^2)
          
        } else {
          
          print(gene)
          print(method)
          
          return(NA)
        }
        
        
      } else {
        
        return(NA)
        
      }
      
      
    }
    
    est_R2 = unlist(sapply(methods, get_R2_one_method))
    
    return(est_R2)
    
  } else {
    
    return(rep(NA, length(methods)))
    
  }
  
  
}

res.mat = t(sapply(genes, get_R2))
colnames(res.mat) = methods

# create empty file to save results
out.dir = file.path(GReX.dir, 
                    paste0("1_Simu_R2_cp_", cp,"_He2_", he2, ".txt"))
out.col = colnames(res.mat)
out.df = data.frame(out.col)
myData = data.frame(matrix(nrow = 0, ncol = length(out.col)))
colnames(myData) = out.col
write.table(myData, 
            out.dir,
            row.names = F,
            col.names = T,
            quote = F,
            sep = '\t')


# write results into table
write.table(res.mat, 
            out.dir,
            row.names = F,
            col.names = F,
            quote = F,
            sep = '\t',
            append = T)



            