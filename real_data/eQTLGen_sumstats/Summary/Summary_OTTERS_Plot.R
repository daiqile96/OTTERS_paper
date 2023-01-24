library(dplyr)
library(ACAT)
res.dir = '~/projects/OTTERS/Data/RDA/Results/'

##################### Read in Predicted GReX for PRS-CS ###### 
GReX = read.table(file.path(res.dir, "GReX_R2.txt"), sep = '\t', header = T)
methods = c("P0.001", "P0.05", "SDPR", "lassosum", "PRScs")
GReX = GReX[rowSums(is.na(GReX[, methods]) == F) > 0,]

compare_R2_tab = function(res, method1, method2){
  
  sub_res = res[, c("TargetID", method1, method2)]
  sub_res = sub_res[complete.cases(sub_res), ]
  colnames(sub_res) = c("TargetID", "M1", "M2")
  maxR = max(sub_res[, c("M1", "M2")])
  
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
    theme(text = element_text(size=10, face = "bold"),
          axis.text = element_text(face="bold")) 
  
  return(p1)
  
}

library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
p1 = compare_R2_tab(GReX, "PRScs", "P0.001") + labs(x = "PRS-CS",y = "P+T(0.001)") + labs(tag = "A")
p2 = compare_R2_tab(GReX, "PRScs", "P0.05") + labs(x = "PRS-CS",y = "P+T(0.05)") + labs(tag = "B")
p3 = compare_R2_tab(GReX, "PRScs", "lassosum") + labs(x = "PRS-CS",y = "lassosum") + labs(tag = "C")
p4 = compare_R2_tab(GReX, "PRScs", "SDPR") + labs(x = "PRS-CS", y = "SDPR")  + labs(tag = "D")
p5 = compare_R2_tab(GReX, "PRScs", "FUSION") + labs(x = "PRS-CS", y = "FUSION")  + labs(tag = "E")


pdf(file.path(res.dir,"Figure3.pdf"), width = 8.8, height = 6)
grid.arrange(p1,p2,p3,p4,p5, nrow=2,
             top = textGrob(bquote("Training"~R^2),gp=gpar(fontsize=20,font=2),
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
  return(t(c(n_test, gt_0.01, med_among_gt_0.01)))
  
}

method_list = c("PRScs", "SDPR", "P0.001", "P0.05", "lassosum", "FUSION")
name_list = c("PRScs", "SDPR", "P+T(0.001)", "P+T(0.05)", "lassosum", "FUSION")
res = sapply(1:length(method_list), function(i) summarize_res(GReX,method_list[i]))
colnames(res) = name_list
rownames(res) = c("Total Genes",
                  "Genes with $R^2$ > 0.01",
                  "Median $R^2$")
pander(res)

##################### eQTLGen ############################
# full.twas.dir = file.path(res.dir, 'TWAS', 'TWAS.txt')
# full.twas = read.table(full.twas.dir, sep = '\t', header = T)
# twas = full.twas[is.na(full.twas$ACAT) == F,]
# the above lines are used to check the number of genes in TWAS


full.twas.dir = file.path(res.dir, 'TWAS', 'TWAS_Filter.txt')
full.twas = read.table(full.twas.dir, sep = '\t', header = T)
twas = full.twas[is.na(full.twas$ACAT) == F,]
bon_P_cut = 0.05/16678
# genome.cutoff = 2.5 * 10^-6

# get significant genes (at least significant in one method)
methods = c('P0.001', 'P0.05', 'lassosum', 'SDPR', 'PRScs')
sig = rowSums(twas[, methods] < bon_P_cut, na.rm = T) > 0
twas.sig.tot.ID = twas$TargetID[sig]
length(twas.sig.tot.ID)
acat.sig = twas[twas[, 'ACAT'] < bon_P_cut, ]

#################### Predicted GReX for all significant genes ########
# GReX_PRScs = read.table(file.path(res.dir, "GReX/PRScs_GReX.txt"), sep = '\t', header = T)
# sub.GReX.cor = cor(t(GReX_PRScs[twas.sig.tot.ID, ]))
# diag(sub.GReX.cor) = 0
# write.table(sub.GReX.cor, file.path(res.dir, "GReX/PRScs_GReX_cor_sig_eQTLGen_TWAS.txt"),
#             col.names = T,
#             row.names = T,
#             sep = '\t',
#             quote = F)

#################### Independent Significant TWAS genes ########
sub.GReX.cor = read.table(file.path(res.dir, "GReX/PRScs_GReX_cor_sig_eQTLGen_TWAS.txt"),
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

full.methods = c("P0.001", "P0.05", "PRScs", "SDPR", "lassosum", "ACAT")
idp.full = lapply(full.methods, get_idp_twas_sig, cutoff = 0.5)
n_idp = lapply(idp.full, nrow) %>% unlist()
n_idp.df = data.frame(Method = full.methods, "Number of Indepdent TWAS Genes" = n_idp, check.names = F)
openxlsx::write.xlsx(n_idp.df,file.path(res.dir, "n_idp.xlsx"))

##################### Manhattan Plot ####################
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
    theme(text = element_text(size = size, face = "bold"), axis.text.x = element_text(angle = -50, face = "bold", vjust=-0.5))
}

get_manhattan_plot = function(TWAS_res, method, size = 40){
  
  
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

png(file.path(res.dir,"Figure4.png"), width = 1400, height = 900, units = "px")
get_manhattan_plot(twas, method = "ACAT", size = 28) + theme(axis.text=element_text(size=20)) + labs(title = "TWAS of Cardiovascular Disease by OTTERS")
dev.off()

pdf(file.path(res.dir,"Figure4.pdf"), width = 18, height = 9)
get_manhattan_plot(twas, method = "ACAT", size = 28) + theme(axis.text=element_text(size=20)) + labs(title = "TWAS of Cardiovascular Disease by OTTERS")
dev.off()

#------------------------- Figure S7 ----------------------------------

library(gridExtra)
png(file.path(res.dir,"FigureS4.png"), width = 2500, height = 1500, units = "px") 
p1 = get_manhattan_plot(twas, "P0.001") + labs(tag = "A") + labs(title = "TWAS of Cardiovascular Disease by P+T(0.001)")
p2 = get_manhattan_plot(twas, "P0.05") + labs(tag = "B") + labs(title = "TWAS of Cardiovascular Disease by P+T(0.05)")
p3 = get_manhattan_plot(twas, "lassosum") + labs(tag = "C") + labs(title = "TWAS of Cardiovascular Disease by lassosum")
p4 = get_manhattan_plot(twas, "SDPR") + labs(tag = "D") + labs(title = "TWAS of Cardiovascular Disease by SDPR")
p5 = get_manhattan_plot(twas, "PRScs") + labs(tag = "E") + labs(title = "TWAS of Cardiovascular Disease by PRS-CS")
grid.arrange(p1, p2, p3, p4, p5, nrow = 3)
dev.off()

png(file.path(res.dir,"FigureS4_extend.png"), width = 2500, height = 1500, units = "px") 
p1 = get_manhattan_plot(twas, "P0.001") + labs(tag = "A") + labs(title = "TWAS of Cardiovascular Disease by P+T(0.001)")
p2 = get_manhattan_plot(twas, "P0.05") + labs(tag = "B") + labs(title = "TWAS of Cardiovascular Disease by P+T(0.05)")
p3 = get_manhattan_plot(twas, "lassosum") + labs(tag = "C") + labs(title = "TWAS of Cardiovascular Disease by lassosum")
p4 = get_manhattan_plot(twas, "SDPR") + labs(tag = "D") + labs(title = "TWAS of Cardiovascular Disease by SDPR")
p5 = get_manhattan_plot(twas, "PRScs") + labs(tag = "E") + labs(title = "TWAS of Cardiovascular Disease by PRS-CS")
p6 = get_manhattan_plot(twas, method = "ACAT") + labs(title = "TWAS of Cardiovascular Disease by OTTERS") + labs(tag = "F") 
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)
dev.off()

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
png(file.path(res.dir,"FigureS5.png"), width = 1500, height = 1000, units = "px") 
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)
dev.off()

#------------------------- Table 2 --------------------------------

ACAT_idp = get_idp_twas_sig(method = "ACAT", cutoff = 0.5) %>% 
  arrange(CHROM, GeneStart)
library(openxlsx)
write.xlsx(ACAT_idp, file.path(res.dir, "Table2.xlsx"))

ACAT_all = twas[twas$ACAT < bon_P_cut, ]
remove.ID = ACAT_all$TargetID[ACAT_all$TargetID %in% ACAT_idp$TargetID == F]
idx = sub.GReX.cor > 0.5
for (ID in remove.ID){
  print(paste0(ID, ":", colnames(sub.GReX.cor)[idx[ID,]]))
}

#------------------------- Table S2 ------------------------------
sig_otherm = NULL
for (method in c("P0.001", "P0.05",  "lassosum", "SDPR", "PRScs")){
  temp_tab = get_idp_twas_sig(method = method, cutoff = 0.5) %>% 
    mutate(method = method) %>% 
    arrange(CHROM, GeneStart) %>% 
    filter(TargetID %in% ACAT_idp$TargetID == F)
  sig_otherm = rbind(sig_otherm, temp_tab)
}

write.xlsx(sig_otherm,
           file.path(res.dir, "TableS2.xlsx"))


#-------------------------- Interested Genes ---------------------

selected.genes = c('CLCN6', 'RP11-378A13.1', 'LINC01093', 'SIDT2', 'FES',
                   'ACE', 'EDN3')
selected.genes.anno = twas[twas$GeneName %in% selected.genes, 1:4]
selected.genes.dir = file.path(res.dir, 'SelectGenes')
dir.create(selected.genes.dir)

write.table(selected.genes.anno,
            file.path(selected.genes.dir, 'SelectGenesAnno.txt'),
            quote = F,
            sep = '\t',
            col.names = T,
            row.names = F)

res = NULL
anno_dir = "~/projects/OTTERS/Data/RDA/GTEx_V8/Anno/"

for (chr in c(1:22)){
  
  anno = read.table(file.path(anno_dir, 
                              paste0('GTEx_CHR', chr, '_GeneAnno.txt')),
                    header = T)
  
  for (batch in (1:ceiling(nrow(anno)/100))){
    
    start = (batch - 1)*100 + 1
    end = min(batch*100, nrow(anno))
    sub_anno = anno[start:end, ]
    
    if (sum(selected.genes.anno$TargetID %in% sub_anno$TargetID)){
      tmp = c(chr, batch)
      res = rbind(res, tmp)
      
    }

  }
  
}

write.table(res, 
            file.path(selected.genes.dir, 'train_group.txt'),
            sep = '\t',
            row.names = F,
            col.names = F,
            quote = F)




library(ggpubr)
selected.genes.anno

get_only_legend <- function(plot) {
  
  # get tabular interpretation of plot
  plot_table <- ggplot_gtable(ggplot_build(plot)) 
  
  #  Mark only legend in plot
  legend_plot <- which(sapply(plot_table$grobs, function(x) x$name) == "guide-box") 
  
  # extract legend
  legend <- plot_table$grobs[[legend_plot]]
  
  # return legend
  return(legend) 
}


plot_weight_5methods = function(row, y_start, y_end){
  
  gene = selected.genes.anno$TargetID[row]
  chr.str = paste0("CHR",selected.genes.anno$CHROM[row])
  
  weight_pos_plot <- function(gene, method,
                              point_alpha = 1,
                              GWAS_sig_level=1e-2,
                              method_name
  ){
    
    e.w = read.table(file.path(selected.genes.dir, chr.str, paste(gene, method, 'e.txt', sep = "_")),
                     header = T) %>% 
      select(CHROM, POS, A1, A2, TargetID, ES, snpID)
    
    g.w = read.table(file.path(selected.genes.dir, chr.str, paste(gene, method, 'g.txt', sep = "_")),
                     header = T) %>% 
      select(CHROM, POS, A1, A2, Z, snpID)
    
    cmp.w = merge(g.w, e.w, by = c("CHROM", "POS", "A1", "A2", "snpID")) %>% 
      arrange(POS) %>% 
      mutate(P = 2 * pnorm(q=abs(Z), lower.tail=F),
             logP = -log10(P),
             GWAS_Zscore = ifelse(Z > 0, 'Positive', 'Negative'),
             GWAS_Zscore = factor(GWAS_Zscore, levels = c('Positive', 'Negative'))) 
    
    p = ggplot(
      data = subset(cmp.w, P >= GWAS_sig_level),
      aes(x=POS,
          y=ES)) +
      geom_point(size=1,
                 alpha=point_alpha) +
      geom_point(
        data=subset(cmp.w, P < GWAS_sig_level),
        aes(
          x=POS, 
          y=ES,
          col=logP,
          fill=logP,
          shape=GWAS_Zscore), 
        size=2, 
        alpha=point_alpha)+
      scale_shape_manual(values = c(24,25)) +
      scale_color_gradient(low="#F0E442",high="#E74C3C") +
      scale_fill_gradient(low="#F0E442",high="#E74C3C") +
      theme(axis.title.x=element_blank(),axis.title.y=element_blank())+  # hide axis labels
      ggtitle(paste0(method_name, ' TWAS p-value:', formatC(twas[twas$TargetID == gene, method]))) +
      geom_hline(
        yintercept=0, 
        linetype='dashed', 
        size=0.3, 
        color="black") +
      labs(color="-log10(GWAS p-value)", fill = "-log10(GWAS p-value)", shape = "GWAS Z Score") +
      theme(plot.title=element_text(size=10,face="bold"),
            legend.title = element_text(color="black",size=10,face="bold"))
    
    return(p)
  }
  
  p1 = weight_pos_plot(gene, method = "P0.001", method_name = "P+T(0.001)") + 
    theme(legend.position = "none") +
    ylim(y_start, y_end)
  p2 = weight_pos_plot(gene, method = "P0.05", method_name = "P+T(0.05)") + 
    theme(legend.position = "none") +
    ylim(y_start, y_end)
  p3 = weight_pos_plot(gene, method = "lassosum", method_name = "lassosum") + 
    theme(legend.position = "none") +
    ylim(y_start, y_end)
  p4 = weight_pos_plot(gene, method = "SDPR", method_name = "SDPR") + 
    theme(legend.position = "none") +
    ylim(y_start, y_end)
  p5 = weight_pos_plot(gene, method = "PRScs", method_name = "PRS-CS") + 
    theme(legend.position = "bottom") +
    ylim(y_start, y_end)
  
  legend = get_only_legend(p5)
  com_plots = gridExtra::grid.arrange(p1, p2, p3, p4, p5 + theme(legend.position = "none"),
                                      nrow = 2)
  p = gridExtra::grid.arrange(com_plots, legend, ncol = 1, heights = c(10, 1))
  
  pdf(file.path(selected.genes.dir, paste0(gene,".pdf")), width = 10, height = 8) 
  p.new = annotate_figure(p, 
                          bottom = text_grob(paste("Position on", chr.str, "(Mb)"),face="bold",size=10),
                          left= text_grob(paste('eQTL Weights for', twas[twas$TargetID == gene, 'GeneName']),
                                          face="bold",size=10,rot=90))
  print(p.new)
  dev.off()
  
}

# for ENSG00000149577 (SIDT2)
plot_weight_5methods(4, -0.7, 0.9)
# for ENSG00000124205 (EDN3)
plot_weight_5methods(7, -0.1, 0.18)
# for ENSG00000249173 (LINC01093) 
plot_weight_5methods(3, -0.22, 0.15)

####################### For FULL TWAS ####################################
full.twas.dir = file.path(res.dir, 'TWAS', 'TWAS.txt')
full.twas = read.table(full.twas.dir, sep = '\t', header = T)
twas = full.twas[is.na(full.twas$ACAT) == F,]
bon_P_cut = 0.05/nrow(GReX)


ACAT_all.full = twas[twas$ACAT < bon_P_cut, ]
sum(ACAT_all.full$TargetID %in% ACAT_all$TargetID == F)

ACAT_all.full[ACAT_all.full$TargetID %in% ACAT_all$TargetID == F, ]
