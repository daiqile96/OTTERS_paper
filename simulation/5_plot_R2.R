.libPaths("/mnt/YangFSS/data/qdai/.local/R/4.0.4/lib64/R/library/")
library(dplyr)
library(ggplot2)
library(ACAT)
Manu_dir = "~/YangFSSdata2/qdai/BS_TWAS/Scripts/Revision"
Simu.dir = "~/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu"

#-------------------------Plot R^2 Box Plot---------------------------

temp_R2_tab = function(res_dir, method, cp, he2){
  
  sub.dir = paste0("std_cp_", cp, "_He2_", he2)
  GReX.dir = file.path(Simu.dir, "Results", sub.dir, "GReX")
  
  r2_matrix_dir = file.path(GReX.dir, 
                            paste0("1_Simu_R2_cp_", cp,"_He2_", he2, ".txt"))
  
  r2_matrix = read.table(r2_matrix_dir, header = T, check.names = F)
  r2_vec = r2_matrix[, method]
  r2_vec[is.na(r2_vec)] = 0
  
  tab = r2_matrix %>% mutate(Method = method, CausalProp = cp, Heritability = he2) %>% 
    mutate(CausalProp = paste0("Causal Proportion=", CausalProp), 
           He2 = factor(Heritability, levels = c(0.01, 0.05, 0.1),
                        labels = c(0.01, 0.05, 0.1)),
           R2 = r2_vec) %>% 
    dplyr::select(CausalProp, He2, Method, R2) 
  
  return(tab)
}

res_dir = "~/YangFSSdata2/qdai/BS_TWAS/Results/ReSimu/"

res_tab = NULL

for (cp in c(0.001, 0.01)){
  for (he2 in c(0.01, 0.05, 0.1)){
    for (method in c("SDPR", "PRScs", "lassosum", "P0.05", "P0.001", "FUSION")){
      temp_tab =  temp_R2_tab(res_dir, method, cp, he2)
      res_tab = rbind(res_tab, temp_tab)
    }
  }
}

res_tab$Method = factor(res_tab$Method,
                        levels = c("P0.001", "P0.05",
                                   "lassosum",  "SDPR", "PRScs", "FUSION"),
                        labels = c("P+T(0.001)", "P+T(0.05)",
                                   "lassosum", "SDPR", "PRS-CS", "FUSION"))

maxR = max(res_tab$R2, na.rm = T)
MyPal = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#89cff0", "#F0027F")

Figure2A = ggplot(res_tab, aes(x = He2, y = R2)) +  # ggplot function
  geom_boxplot(aes(fill=Method),
               position=position_dodge(.7), width = 0.5,
               outlier.size = 0.5) + 
  facet_wrap(~CausalProp, ncol = 2) +
  scale_y_continuous(trans = "sqrt", breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.5), limits = c(0, maxR)) +
  theme(text = element_text(size=18,face="bold"),
        legend.position="bottom",
        legend.box = "vertical",
        axis.text = element_text(face="bold")) +
  labs(y = bquote("Test"~R^2), x = bquote(h[e]^2))  +
  scale_fill_manual(values=MyPal)

pdf(file.path(Manu_dir,"Diff_Genes_Simu_R2.pdf"), width = 7, height = 6)
Figure2A + theme( legend.position = "none" )
dev.off()


############### POWER #################

power.dir = file.path(Simu.dir, "Results", "2_ReSimu_Power_100_per_group.txt")
power.tab = data.frame(read.table(power.dir, header = F))
colnames(power.tab) = c("CausalProp", "He2", "N", "SDPR", "PRScs", "lassosum", "P0.05", "P0.001", "FUSION", "ACAT")

idx = power.tab$He2 == 0.01 & power.tab$N %in% c(50, 75, 100, 150)
power.tab = power.tab[!idx, ]

temp_tab = function(method){
  
  power.tab$Power = power.tab[, method]
  
  tab = power.tab %>% mutate(Method = method) %>% 
    mutate(CausalProp = paste0("Causal Proportion=", CausalProp), 
           He2 = paste0("Expression Heritabillity=", He2)) %>% 
    dplyr::select(CausalProp, Method, Power, He2, N) 
  
  return(tab)
}

power_df_plot <- rbind(temp_tab("SDPR"),
                       temp_tab("PRScs"),
                       temp_tab("P0.05"),
                       temp_tab("P0.001"),
                       temp_tab("lassosum"),
                       temp_tab("FUSION"),
                       temp_tab("ACAT"))

power_df_plot$Method = factor(power_df_plot$Method,
                              levels = c("P0.001", "P0.05",
                                         "lassosum", "SDPR", "PRScs", "FUSION", "ACAT"),
                              labels = c("P+T(0.001)", "P+T(0.05)",
                                         "lassosum",  "SDPR", "PRS-CS","FUSION", "OTTERS"))

MyPal = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#89cff0", "#F0027F")

Figure2B = ggplot(power_df_plot, aes(x = factor(N), y = as.numeric(Power),
                                     fill = Method)) +
  geom_col(width = 0.6, position = position_dodge(0.7)) +
  facet_wrap(~CausalProp + He2, ncol = 3, scales = "free_x") +
  labs(y = "TWAS Power", x = "n_gwas (in K)", title = "") +
  theme(text = element_text(size=18, face="bold"),
        legend.position="bottom",
        legend.box = "vertical",
        axis.text = element_text(face="bold")) +
  scale_fill_manual(values=MyPal)

pdf(file.path(Manu_dir,"Diff_Genes_Simu_Power.pdf"), width = 11, height = 6)
Figure2B + theme( legend.position = "none" )
dev.off()


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

legend = get_only_legend(Figure2B)

png(file.path(Manu_dir,"Figure2_legend.pdf"), width = 18, height = 1, units='mm', res = 300)
plot(legend)
dev.off()


com_plots = gridExtra::grid.arrange(Figure2A + theme( legend.position = "none" ) + labs(tag = "A"), 
                                    Figure2B + theme( legend.position = "none" ) + labs(tag = "B"), 
                                    widths = c(4.4, 5.6),
                                    nrow = 1)

pdf(file.path(Manu_dir,"Figure2.pdf"), width = 17, height = 6)
gridExtra::grid.arrange(com_plots, legend, ncol = 1, heights = c(10, 1))
dev.off()

png(file.path(Manu_dir,"Figure2.png"), width = 700, height = 300, units='mm', res = 300)
gridExtra::grid.arrange(com_plots, legend, ncol = 1, heights = c(10, 1))
dev.off()


############### POWER_by_group #################

power.dir = file.path(Simu.dir, "Results", "2_ReSimu_Power_by_group_100_per_group.txt")
power.tab = data.frame(read.table(power.dir, header = F))
colnames(power.tab) = c("CausalProp", "He2", "N", "Group", "SDPR", "PRScs", "lassosum", "P0.05", "P0.001", "FUSION", "ACAT")

idx = power.tab$He2 == 0.01 & power.tab$N %in% c(50, 75, 100, 150)
power.tab = power.tab[!idx, ]

temp_tab = function(method){
  
  power.tab$Power = power.tab[, method]
  
  tab = power.tab %>% mutate(Method = method) %>% 
    mutate(CausalProp = paste0("Causal Proportion=", CausalProp), 
           He2 = factor(He2, levels = c(0.01,0.05,0.1),
                        labels = c(0.01,0.05,0.1))) %>% 
    dplyr::select(CausalProp, Method, Group, Power, He2, N) 
  
  return(tab)
}

power_df_plot <- rbind(temp_tab("SDPR"),
                       temp_tab("PRScs"),
                       temp_tab("P0.05"),
                       temp_tab("P0.001"),
                       temp_tab("lassosum"),
                       temp_tab("FUSION"),
                       temp_tab("ACAT"))

power_df_plot$Method = factor(power_df_plot$Method,
                              levels = c("P0.001", "P0.05",
                                         "lassosum", "SDPR", "PRScs", "FUSION", "ACAT"),
                              labels = c("P+T(0.001)", "P+T(0.05)",
                                         "lassosum",  "SDPR", "PRS-CS","FUSION", "OTTERS"))

MyPal = c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#89cff0", "#F0027F")

sum_df = power_df_plot %>% group_by(CausalProp, He2, N, Method) %>% summarise(meanPower = mean(Power))

plot_by_senario = function(he2){
  
  sub_df = power_df_plot %>% filter(He2 == he2) %>% 
    mutate(Group = paste("Group", Group)) 
  ggplot(sub_df, aes(x = factor(N), y = as.numeric(Power),
                     fill = Method)) +
    geom_col(width = 0.6, position = position_dodge(0.7)) +
    facet_wrap(~CausalProp+Group, ncol = 5) +
    labs(y = "TWAS Power", x = "n_gwas (in K)", title = "") +
    theme(text = element_text(size=20, face="bold"),
          legend.position="bottom",
          legend.box = "vertical",
          axis.text = element_text(face="bold")) +
    scale_fill_manual(values=MyPal)
  
}


pdf(file.path(Manu_dir,"Diff_Genes_Simu_Power_he20.01.pdf"), width = 21, height = 9)
plot_by_senario(0.01)
dev.off()

pdf(file.path(Manu_dir,"Diff_Genes_Simu_Power_he20.05.pdf"), width = 21, height = 9)
plot_by_senario(0.05)
dev.off()


pdf(file.path(Manu_dir,"Diff_Genes_Simu_Power_he20.1.pdf"), width = 21, height = 9)
plot_by_senario(0.1)
dev.off()


