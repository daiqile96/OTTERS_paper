# Settings
require(dplyr, quietly = T, warn.conflicts = F)
set.seed(2022)

# Read in gene annotation file (in hg19)
wk.dir = "~/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu"
GeneLength = read.table(file.path(wk.dir, "0_GeneLength.txt"), header = T)
GeneLength = GeneLength %>% 
  filter(CHROM %in% c(1:22)) %>% 
  mutate(Length = GeneEnd - GeneStart) %>% 
  arrange(Length)

# select by quantile
length.qtl = quantile(GeneLength$Length, probs = seq(0, 1, 0.2))

# select 20 by group 
# sample.by.length = function(qtl){
#   
#   sub.GeneLength = GeneLength %>% 
#     filter(Length < length.qtl[qtl] & Length > length.qtl[qtl - 1]) 
#   
#   select.idx = sample(1:nrow(sub.GeneLength), 20)
#   
#   return(sub.GeneLength[select.idx, 'TargetID'])
#   
# }
# 
# # return the list of selected genes
# select.genes = unlist(lapply(2:6, sample.by.length))
# 
# # extract the annotation table for the selected genes
# res.tab = GeneLength %>% filter(TargetID %in% select.genes)
# 
# # write the table to file
# write.table(res.tab, file = file.path(wk.dir, "0_ReSimu_GeneAnno.txt"), 
#             row.names = F,
#             col.names = F,
#             quote = F)
# 
# print("Done Selecting Genes.")


sample.by.length = function(qtl){
  
  sub.GeneLength = GeneLength %>% 
    filter(Length < length.qtl[qtl] & Length > length.qtl[qtl - 1]) 
  
  select.idx = sample(1:nrow(sub.GeneLength), 100)
  
  return(sub.GeneLength[select.idx, 'TargetID'])
  
}

# return the list of selected genes
select.genes = unlist(lapply(2:6, sample.by.length))

# extract the annotation table for the selected genes
res.tab = GeneLength %>% filter(TargetID %in% select.genes)

# write the table to file
write.table(res.tab, file = file.path(wk.dir, "0_ReSimu_GeneAnno_100_per_group.txt"), 
            row.names = F,
            col.names = F,
            quote = F)

print("Done Selecting Genes.")
