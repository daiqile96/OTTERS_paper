################ paramters that can be changed ############
methods = c('P0.001', 'P0.05', 'lassosum', 'SDPR', 'PRScs', 'FUSION')

# directory of results
sub.res.dir = 'ReRDA/Results'

################## True Expression ##################
data.dir = '~/projects/OTTERS/Data'
exp.dir = file.path(data.dir, "RDA/GTEx_V8", 'Whole_Blood_GTEx_Exp.txt')
Exp = read.table(exp.dir, sep = '\t', header = T, check.names = F)
Exp.mat = as.matrix(Exp[, 6:ncol(Exp)])
Gene.ID = Exp$TargetID
Sub.ID = colnames(Exp.mat)
rownames(Exp.mat) = Gene.ID
print('Done Reading True Expression')

################## Job List #######################
# Batch list 
jobs = read.table(file.path(data.dir, 'RDA/GTEx_V8/Anno/train_group.txt'), sep = '\t',
                  header = F)
colnames(jobs) = c('chr', 'batch')

################## Predicted GReX #######################
GReX_Method = function(method){
  
  method.str = paste0(method, '.txt')
  
  test.R2 = NULL
  pred.GReX.full = NULL
  
  for (row in 1:nrow(jobs)){
    
    batch = jobs[row, 'batch']
    chr = jobs[row, 'chr']
    
    batch.str = paste0('Batch', batch)
    chr.str = paste0('CHR', chr)
    
    Pred.GReX.dir = file.path(data.dir, sub.res.dir, chr.str, batch.str, 'GReX', method.str)
    Pred.GReX = read.table(Pred.GReX.dir, header = T, sep = '\t',
                           skip = 1, row.names = 1, check.names = F)
    
    inter.genes = intersect(row.names(Pred.GReX), Gene.ID)
    Pred.GReX = Pred.GReX[inter.genes, Sub.ID]
    Sub.Exp.mat = Exp.mat[inter.genes, Sub.ID]
    
    test.R2.tmp = sapply(seq.int(dim(Pred.GReX)[1]), function(i) cor(unlist(Pred.GReX[i,]), Sub.Exp.mat[i,])^2)
    test.R2.tmp.df = data.frame(inter.genes, test.R2.tmp)
    
    test.R2 = rbind(test.R2, test.R2.tmp.df)
    pred.GReX.full = rbind(pred.GReX.full, Pred.GReX)
  }
  
  colnames(test.R2) = c('TargetID', method)
  out.dir = file.path(data.dir, sub.res.dir, 'GReX')
  if (dir.exists(out.dir) == F){
    dir.create(out.dir, recursive = T)
  }
  
  write.table(test.R2, file.path(out.dir, method.str),
              quote = F, col.names = T, row.names = F,
              sep = '\t')
  
  if (method == 'lassosum'){
    write.table(pred.GReX.full, file.path(out.dir, paste0(method, '_GReX.txt')),
                quote = F, col.names = T, row.names = T,
                sep = '\t')
  }
  
  print(paste('Finish', method))
  print(paste('Out file:', file.path(out.dir, method.str)))
  
}

# Get test R2 for each method
print('Start calculating test R2')
sapply(methods, GReX_Method)
print('Done')

# Get test R2 for all methods
print('Start summarizing test R2')
GReX = NULL
for (method in methods){
  
  method.str = paste0(method, '.txt')
  test.R2 = read.table(file.path(data.dir, sub.res.dir, 'GReX', method.str),
                       sep = '\t',
                       header = T)
  test.R2 = data.frame(test.R2)
  
  if (is.null(GReX)){
    
    GReX = test.R2
    
  } else {
    
    GReX = merge(GReX, test.R2, all = T)
    
  }
  
}

write.table(GReX, file.path(data.dir, sub.res.dir, 'GReX_R2.txt'),
            sep = '\t',
            col.names = T,
            quote = F)
print('Done')



