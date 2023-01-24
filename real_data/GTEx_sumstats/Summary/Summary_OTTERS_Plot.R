library(dplyr)
library(ACAT)
res.dir = '~/projects/OTTERS/Data/ReRDA/Results/'
methods = c("lassosum", "P0.001", "P0.05", "PRScs", "SDPR", "ACAT")

##################### Read in Predicted GReX for PRS-CS ###### 
GReX = read.table(file.path(res.dir, "GReX_R2.txt"), sep = '\t', header = T)

# compare_R2_tab = function(res, method1, method2){
#   
#   sub_res = res[, c("TargetID", method1, method2)]
#   sub_res = sub_res[complete.cases(sub_res), ]
#   colnames(sub_res) = c("TargetID", "M1", "M2")
#   maxR = max(sub_res[, c("M1", "M2")])
#   
#   sub_res = sub_res %>% mutate(col = case_when(
#     M1 >=0.01 & M2 <0.01 ~ "g1",
#     M1 <0.01 & M2 >=0.01 ~ "g2",
#     M1 >=0.01 & M2 >=0.01 ~ "g3",
#   )) %>% filter(M1 >=0.01 | M2 >= 0.01)
#   
#   cols = c("#F8766D", "#619CFF")
#   p1 = ggplot(sub_res, aes(x=M1, y=M2, color = col)) +
#     geom_point(show.legend = FALSE, size = 1) +
#     scale_x_continuous(trans = "sqrt", breaks = c(0,0.01,0.1,0.2, 0.5), limits = c(0, maxR)) +
#     scale_y_continuous(trans = "sqrt", breaks = c(0,0.01,0.1,0.2, 0.5), limits = c(0, maxR)) +
#     geom_abline(intercept = 0, slope = 1, color="black", linetype = "dotdash") +
#     labs(x = method1, y = method2) +
#     theme(text = element_text(size=20, face = "bold"),
#           axis.text = element_text(face="bold")) +
#     scale_color_manual(values=cols[(3-length(unique(sub_res$col))):2])
#   
#   return(p1)
#   
# }

compare_R2_tab = function(res, method1, method2){

  sub_res = res[, c("TargetID", method1, method2)]
  sub_res = sub_res[complete.cases(sub_res), ]
  colnames(sub_res) = c("TargetID", "M1", "M2")
  maxR = max(sub_res[, c("M1", "M2")])

  cols <- c("g1" = "#F8766D", "g2" = "#00BA38", "g3" = "#619CFF")
  sub_res = sub_res %>% mutate(col = case_when(
    M1 >=0.01 & M2 <0.01 ~ "g1",
    M1 <0.01 & M2 >=0.01 ~ "g2",
    M1 >=0.01 & M2 >=0.01 ~ "g3",
  )) %>% filter(M1 >=0.01 | M2 >= 0.01)
  
  p1 = ggplot(sub_res, aes(x=M1, y=M2, color = col)) +
    geom_point(show.legend = FALSE, size = 1) +
    scale_x_continuous(trans = "sqrt", breaks = c(0,0.01,0.1,0.2, 0.5), limits = c(0, maxR)) +
    scale_y_continuous(trans = "sqrt", breaks = c(0,0.01,0.1,0.2, 0.5), limits = c(0, maxR)) +
    geom_abline(intercept = 0, slope = 1, color="black", linetype = "dotdash") +
    labs(x = method1, y = method2) +
    theme(text = element_text(size=20, face = "bold"),
          axis.text = element_text(face="bold")) +
    scale_color_manual(values = cols)

  return(p1)

}

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
# p1 = compare_R2_tab(GReX, "PRScs", "P0.001") + labs(x = "PRS-CS",y = "P+T(0.001)") + labs(tag = "A")
# p2 = compare_R2_tab(GReX, "PRScs", "P0.05") + labs(x = "PRS-CS",y = "P+T(0.05)") + labs(tag = "B")
# p3 = compare_R2_tab(GReX, "PRScs", "lassosum") + labs(x = "PRS-CS",y = "lassosum") + labs(tag = "C")
# p4 = compare_R2_tab(GReX, "PRScs", "SDPR") + labs(x = "PRS-CS", y = "SDPR")  + labs(tag = "D")
# p5 = compare_R2_tab(GReX, "PRScs", "FUSION") + labs(x = "PRS-CS", y = "FUSION")  + labs(tag = "E")


p1 = compare_R2_tab(GReX, "lassosum", "P0.001") + labs(x = "lassosum",y = "P+T(0.001)") + labs(tag = "A")
p2 = compare_R2_tab(GReX, "lassosum", "P0.05") + labs(x = "lassosum",y = "P+T(0.05)") + labs(tag = "B")
p3 = compare_R2_tab(GReX, "lassosum", "PRScs") + labs(x = "lassosum",y = "PRS-CS") + labs(tag = "C")
p4 = compare_R2_tab(GReX, "lassosum", "SDPR") + labs(x = "lassosum", y = "SDPR")  + labs(tag = "D")
p5 = compare_R2_tab(GReX, "lassosum", "FUSION") + labs(x = "lassosum", y = "FUSION")  + labs(tag = "E")


pdf(file.path(res.dir,"Figure3.pdf"), width = 14, height = 9)
grid.arrange(p1,p2,p3,p4,p5, nrow=2,
             top = textGrob(bquote("Training"~R^2),gp=gpar(fontsize=32,font=2),
                            x = 0.08, hjust = 0))
dev.off()

png(file.path(res.dir,"Figure3.png"), width = 140*2.5, height = 90*2.5, units = "mm",
    res = 300)
grid.arrange(p1,p2,p3,p4,p5, nrow=2,
             top = textGrob(bquote("Training"~R^2),gp=gpar(fontsize=32,font=2),
                            x = 0.08, hjust = 0))
dev.off()

library(pander)
summarize_res = function(res_all, method){
  
  res_vec = unlist(res_all[, method])
  res_vec = as.numeric(res_vec[is.na(res_vec) == FALSE])
  
  n_test = length(res_vec)
  gt_0.01 = as.integer(sum(res_vec > 0.01))
  med_among_gt_0.01 = median(res_vec[res_vec > 0.01])
  # gt_0.05 = as.integer(sum(res_vec > 0.05))
  # med_among_gt_0.05 = median(res_vec[res_vec > 0.05])
  return(t(c(n_test, gt_0.01, med_among_gt_0.01, 19696 - gt_0.01)))
  
}

method_list = c("PRScs", "SDPR", "P0.001", "P0.05", "lassosum", "FUSION")
name_list = c("PRScs", "SDPR", "P+T(0.001)", "P+T(0.05)", "lassosum", "FUSION")
res = sapply(1:length(method_list), function(i) summarize_res(GReX,method_list[i]))
colnames(res) = name_list
rownames(res) = c("Total Genes",
                  "Genes with $R^2$ > 0.01",
                  "Median $R^2$",
                  "Missed Genes")
pander(res)


##################### eQTLGen ############################
full.twas.dir = file.path(res.dir, 'TWAS', 'TWAS_Filter.txt')
full.twas = read.table(full.twas.dir, sep = '\t', header = T)
twas = full.twas[is.na(full.twas$ACAT) == F,]
bon_P_cut = 0.05/nrow(GReX)
# genome.cutoff = 2.5 * 10^-6
  
# get significant genes (at least significant in one method)
methods = c('P0.001', 'P0.05', 'lassosum', 'SDPR', 'PRScs', 'FUSION')
sig_tab = twas[, c(methods, 'ACAT')] < bon_P_cut
all_sig = colSums(sig_tab, na.rm = T)
cmp_fusion_acat = sum(sig_tab[, 'FUSION'] == T & sig_tab[, 'ACAT'] == T)
sig = rowSums(twas[, methods] < bon_P_cut, na.rm = T) > 0
twas.sig.tot.ID = twas$TargetID[sig]
acat.sig = twas[twas[, 'ACAT'] < bon_P_cut, ]

#################### Predicted GReX for all significant genes ########
# GReX_lassosum = read.table(file.path(res.dir, "GReX/lassosum_GReX.txt"), sep = '\t', header = T)
# sub.GReX.cor = cor(t(GReX_lassosum[twas.sig.tot.ID, ]))
# diag(sub.GReX.cor) = 0
# write.table(sub.GReX.cor, file.path(res.dir, "GReX/lassosum_GReX_cor_sig_eQTLGen_TWAS.txt"),
#             col.names = T,
#             row.names = T,
#             sep = '\t',
#             quote = F)

#################### Independent Significant TWAS genes ########
sub.GReX.cor = read.table(file.path(res.dir, "GReX/lassosum_GReX_cor_sig_eQTLGen_TWAS.txt"),
                          header = T,
                          row.names = 1,
                          sep = '\t')

get_idp_twas_sig = function(method, cutoff){
  twas = twas[complete.cases(twas[, method]), ]
  twas.sig.method = twas[twas[, method] < bon_P_cut, ] 
  order.sig = order(twas.sig.method[, method])
  twas.sig.method = twas.sig.method[order.sig, ]
  ID.order.sig = unlist(twas.sig.method$TargetID)
  sub.GReX.cor = sub.GReX.cor[ID.order.sig, ID.order.sig]
  
  ID.keep = 1
  for (row in 2:nrow(twas.sig.method)){
    
    if (sum(sub.GReX.cor[row, 1:row]^2 >  cutoff) == 0){
      ID.keep = c(ID.keep, row)
    } 
  }
  
  twas.sig.method[ID.keep, ]
}

full.methods = c("P0.001", "P0.05", "PRScs", "SDPR", "lassosum", "FUSION", "ACAT")
idp.full = lapply(full.methods, get_idp_twas_sig, cutoff = 0.5)
n_idp = lapply(idp.full, nrow) %>% unlist()
n_idp.df = data.frame(Method = full.methods, "Number of Indepdent TWAS Genes" = n_idp, check.names = F)
openxlsx::write.xlsx(n_idp.df,file.path(res.dir, "n_idp.xlsx"))

##################### Manhattan Plot ####################
methods = c("P0.001", "P0.05", "PRScs", "SDPR", "lassosum", "FUSION")
max_p = max(-log10(twas[, c(methods, "ACAT")]), na.rm = T)

myManPlot <- function(manPlot_dt, title = "Manhantton Plot", chr_vec = 1:22, ntop = 3, sig_level = 2.5e-6, size = 23, chrGAP = 5e2, label_genes, max_p){
  # manPlot is a data frame containing columns: CHR, POS, PVALUE, ID
  # CHR column should be of the format: chr*
  # Setup ploting positions
  # chr_vec = sort(unique(manPlot$CHR))
  endPos = 0;
  plotPos = NULL; temp_dt = NULL;
  chrEnd = NULL; LabBreaks = NULL; xlabels = NULL;
  for (chr in chr_vec) {
    print(chr)
    temp = manPlot_dt[manPlot_dt$CHR == chr, ]
    temp$POS = order(temp$POS)
    if(nrow(temp)>0){
      temp_dt = rbind(temp_dt, temp)
      chrPos = (temp$POS - min(temp$POS, na.rm = TRUE) ) + endPos + 1
      endPos = max(chrPos, na.rm = TRUE) + chrGAP
      plotPos = c(plotPos, chrPos)
      yline_pos = max(chrPos, na.rm = TRUE) + chrGAP/2
      chrEnd = c(chrEnd, yline_pos)
      LabBreaks = c(LabBreaks, mean(chrPos, na.rm = TRUE) )
      xlabels = c(xlabels, chr)
    }
  }
  
  manPlot_dt_sort = data.frame(plotPos = plotPos, temp_dt)
  
  require("ggrepel")
  ggplot(manPlot_dt_sort, aes(plotPos, -log10(PVALUE), color = factor(CHR) )) +
    geom_point() + guides(color = FALSE) +
    ylim(c(0, max(-log10(sig_level), max_p))) +
    geom_hline(yintercept = -log10(sig_level), color = "red") +
    labs(x = "Chromosome", y = "-log10(PVALUE)", title = title) +
    scale_x_continuous(breaks = LabBreaks, labels = xlabels, limits = range(plotPos)) +
    geom_label_repel(data = manPlot_dt_sort[manPlot_dt_sort$ID %in% label_genes, ],
                     aes(plotPos, -log10(PVALUE), label = ID), size = 5,
                     nudge_x = .15,
                     box.padding = unit(0.25, "lines"),
                     nudge_y = 1,
                     segment.curvature = -0.1,
                     segment.ncp = 3,
                     segment.angle = 20,
                     max.overlaps =20) +
    theme(text = element_text(size = size, face = "bold"), 
          axis.text.x = element_text(angle = -50, face = "bold", vjust=-0.5))
}


get_manhattan_plot = function(TWAS_res, method, size = 20){
  
  
  sub_TWAS = TWAS_res[, c("CHROM", "GeneStart", method, "GeneName")]
  sub_TWAS = sub_TWAS[complete.cases(sub_TWAS), ] 
  colnames(sub_TWAS) = c("CHR", "POS", "PVALUE", "ID") 
  
  sub_TWAS = sub_TWAS %>% mutate(
    
    CHR = as.numeric(CHR),
    POS = as.numeric(POS)
  ) %>% arrange(CHR, POS)
  
  # sig_level = 0.05/nrow(sub_TWAS)
  genes = get_idp_twas_sig(method, 0.5)$GeneName
  
  print(genes)
  
  p1 = myManPlot(sub_TWAS, 
                 title = method, 
                 chr_vec = c(1:22),
                 size = size, chrGAP = 1.6e2,
                 label_genes = genes,
                 max_p = max_p)
  
  return(p1)
}

# png(file.path(res.dir,"OTTERS_TWAS.png"), width = 1100, height = 600)
# get_manhattan_plot(twas, method = "ACAT", size = 25) + theme(axis.text=element_text(size=18)) + labs(title = "TWAS of Cardiovascular Disease by OTTERS")
# dev.off()

FigureA = get_manhattan_plot(twas, method = "ACAT", size = 25) + 
  theme(axis.text=element_text(size=18)) + 
  labs(title = "TWAS of Cardiovascular Disease by OTTERS")

FigureB = get_manhattan_plot(twas, method = "FUSION", size = 25) + 
  theme(axis.text=element_text(size=18)) + 
  labs(title = "TWAS of Cardiovascular Disease by FUSION")

png(file.path(res.dir,"Figure4.png"), 
    width = 140*2.5, height = 180*2.5, units = "mm", res = 400)
grid.arrange(FigureA + labs(tag = "A"), FigureB + labs(tag = "B"),
             nrow = 2)
dev.off()

pdf(file.path(res.dir,"Figure4.pdf"), width = 14, height = 9)
get_manhattan_plot(twas, method = "ACAT", size = 25) + theme(axis.text=element_text(size=18)) + labs(title = "TWAS of Cardiovascular Disease by OTTERS")
dev.off()

#------------------------- Figure S7 ----------------------------------

library(gridExtra)
png(file.path(res.dir,"FigureS4.png"), width = 2500, height = 1500, units = "px") 
p1 = get_manhattan_plot(twas, "P0.001") + labs(tag = "A") + labs(title = "TWAS of Cardiovascular Disease by P+T(0.001)")
p2 = get_manhattan_plot(twas, "P0.05") + labs(tag = "B") + labs(title = "TWAS of Cardiovascular Disease by P+T(0.05)")
p3 = get_manhattan_plot(twas, "lassosum") + labs(tag = "C") + labs(title = "TWAS of Cardiovascular Disease by lassosum")
p4 = get_manhattan_plot(twas, "SDPR") + labs(tag = "D") + labs(title = "TWAS of Cardiovascular Disease by SDPR")
p5 = get_manhattan_plot(twas, "PRScs") + labs(tag = "E") + labs(title = "TWAS of Cardiovascular Disease by PRS-CS")
p6 = get_manhattan_plot(twas, "FUSION") + labs(tag = "E") + labs(title = "TWAS of Cardiovascular Disease by FUSION")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)
dev.off()


# png(file.path(res.dir,"FigureS4_extend.png"), width = 2500, height = 1500, units = "px") 
# p1 = get_manhattan_plot(twas, "P0.001") + labs(tag = "A") + labs(title = "TWAS of Cardiovascular Disease by P+T(0.001)")
# p2 = get_manhattan_plot(twas, "P0.05") + labs(tag = "B") + labs(title = "TWAS of Cardiovascular Disease by P+T(0.05)")
# p3 = get_manhattan_plot(twas, "lassosum") + labs(tag = "C") + labs(title = "TWAS of Cardiovascular Disease by lassosum")
# p4 = get_manhattan_plot(twas, "SDPR") + labs(tag = "D") + labs(title = "TWAS of Cardiovascular Disease by SDPR")
# p5 = get_manhattan_plot(twas, "PRScs") + labs(tag = "E") + labs(title = "TWAS of Cardiovascular Disease by PRS-CS")
# p6 = get_manhattan_plot(twas, method = "ACAT") + labs(title = "TWAS of Cardiovascular Disease by OTTERS") + labs(tag = "F") 
# grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)
# dev.off()

#------------------------------ Figure S8 ----------------------------------

myqq <- function(pvector, title="", size = 24) {
  
  ci <- 0.95
  pvector = pvector[!is.na(pvector)]
  n <- length(pvector)
  print(n)
  plotdata <- data.frame(
    observed = -log10(sort(pvector)),
    expected = -log10(1:n/n),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(n), shape2 = rev(seq(n)))),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(n), shape2 = rev(seq(n))))
  )
  
  ggplot(plotdata, aes(x = expected, y = observed)) +
    geom_point() + 
    geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
    xlim(c(0, max(plotdata$expected))) + 
    ylim(c(0, max_p)) +
    labs(x = "Expected", y = "Observed", title = title) +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    theme(text = element_text(size = size, face = "bold"))
  
}

get_qq_plot = function(TWAS, method){
  
  sub_TWAS = TWAS[, c("CHROM", "GeneStart", method, "GeneName")]
  sub_TWAS = sub_TWAS[complete.cases(sub_TWAS), ]
  colnames(sub_TWAS) = c("CHR", "POS", "PVALUE", "ID")
  
  p2 = myqq(sub_TWAS$PVALUE, title = "")
  
  return(p2)
  
}

p1 = get_qq_plot(twas, "P0.001") + labs(tag = "A") + labs(title = "P+T(0.001)")
p2 = get_qq_plot(twas, "P0.05") + labs(tag = "B") + labs(title = "P+T(0.05)")
p3 = get_qq_plot(twas, "lassosum") + labs(tag = "C") + labs(title = "lassosum")
p4 = get_qq_plot(twas, "SDPR") + labs(tag = "D") + labs(title = "SDPR")
p5 = get_qq_plot(twas, "PRScs") + labs(tag = "E") + labs(title = "PRS-CS")
p6 = get_qq_plot(twas, "ACAT") + labs(tag = "F") + labs(title = "OTTERS")
p7 = get_qq_plot(twas, "FUSION") + labs(tag = "F") + labs(title = "FUSION")
png(file.path(res.dir,"FigureS5.png"), width = 1500, height = 1500, units = "px") 
grid.arrange(p1, p2, p3, p4, p5, p6, p7, nrow = 3)
dev.off()



#------------------------- Table 2 --------------------------------

ACAT_idp = get_idp_twas_sig(method = "ACAT", cutoff = 0.5) %>% 
  arrange(CHROM, GeneStart)
library(openxlsx)
write.xlsx(ACAT_idp, file.path(res.dir, "Table2.xlsx"))

# pander(ACAT_idp[, c("CHROM", "GeneName", methods, "FUSION", "ACAT"), ])

ACAT_all = twas[twas$ACAT < bon_P_cut, ]
remove.ID = ACAT_all$TargetID[ACAT_all$TargetID %in% ACAT_idp$TargetID == F]
idx = sub.GReX.cor > 0.5
for (ID in remove.ID){
  print(paste0(ID, ":", colnames(sub.GReX.cor)[idx[ID,]]))
}

#------------------------- Table S2 ------------------------------
sig_otherm = NULL
for (method in c("P0.001", "P0.05",  "lassosum", "SDPR", "PRScs", "FUSION")){
  temp_tab = get_idp_twas_sig(method = method, cutoff = 0.5) %>% 
    mutate(method = method) %>% 
    arrange(CHROM, GeneStart) %>% 
    filter(TargetID %in% ACAT_idp$TargetID == F)
  sig_otherm = rbind(sig_otherm, temp_tab)
}

write.xlsx(sig_otherm,
           file.path(res.dir, "TableS2.xlsx"))


####################### For FULL TWAS ####################################
full.twas.dir = file.path(res.dir, 'TWAS', 'TWAS.txt')
full.twas = read.table(full.twas.dir, sep = '\t', header = T)
twas = full.twas[is.na(full.twas$ACAT) == F,]
bon_P_cut = 0.05/nrow(GReX)


ACAT_all.full = twas[twas$ACAT < bon_P_cut, ]
sum(ACAT_all.full$TargetID %in% ACAT_all$TargetID == F)

ACAT_all.full[ACAT_all.full$TargetID %in% ACAT_all$TargetID == F, ]


