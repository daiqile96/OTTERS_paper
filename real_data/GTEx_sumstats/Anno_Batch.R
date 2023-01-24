
anno_dir = "~/projects/OTTERS/Data/RDA/GTEx_V8/Anno/"

res = NULL

for (chr in c(1:22)){
  
  anno = read.table(file.path(anno_dir, 
                              paste0('GTEx_CHR', chr, '_GeneAnno.txt')),
                    header = T)
  
  for (batch in (1:ceiling(nrow(anno)/100))){
    
    start = (batch - 1)*100 + 1
    end = min(batch*100, nrow(anno))
    sub_anno = anno[start:end, ]
    
    write.table(sub_anno,
                file = file.path(anno_dir,
                               paste0('Batch', batch,
                                      '_GTEx_CHR', chr, 
                                      '_GeneAnno.txt')),
                sep = '\t',
                row.names = F,
                col.names = T,
                quote = F)
    
    tmp = c(chr, batch)
    res = rbind(res, tmp)
    
  }

}

write.table(res, 
            file.path(anno_dir, 'train_group.txt'),
            sep = '\t',
            row.names = F,
            col.names = F,
            quote = F)
