library(ggplot2)
# read in all the p-values
Simu.dir = "~/YangFSSdata2/qdai/BS_TWAS/Data/Simulation/ReSimu"
Perm.dir = file.path(Simu.dir, '/Results/std_cp_0.001_He2_0.1/TWAS/Perm')
Perm.res = readRDS(file.path(Perm.dir, 'Perm_res.rds'))

ectract_method = function(mat, method) return(mat[, method])
ectract_method_pvector = function(method) return(as.vector(sapply(Perm.res, ectract_method, method = method)))
methods = c("SDPR", "PRScs", "lassosum", "P0.05", "P0.001", "ACAT.res")
p_matrix = sapply(methods, ectract_method_pvector)
p_max = max(-log10(p_matrix), na.rm = T)

myqq <- function(pvector, title="Quantile-Quantile plot of p-values", size = 24) {
  
  
  ci <- 0.95
  pvector = pvector[!is.na(pvector)]
  n <- length(pvector)
  
  plotdata <- data.frame(
    observed = -log10(sort(pvector)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(n), shape2 = rev(seq(n)))),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(n), shape2 = rev(seq(n))))
  )
  
  ggplot(plotdata, aes(x = expected, y = observed)) +
    geom_point() + 
    geom_ribbon(aes(x = expected, ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
    labs(x = "Expected", y = "Observed", title = title) +
    geom_abline(intercept = 0, slope = 1, colour = "red") +
    theme(text = element_text(size = size, face = "bold")) +
    coord_cartesian(ylim = c(0, p_max), xlim = c(0, max(plotdata$expected))) 
  
}

library(gridExtra)
out.dir = '~/YangFSSdata2/qdai/BS_TWAS/Scripts/Revision'
png(file.path(out.dir, "ReSimu_Perm_QQ.png"), width = 1000, height = 800, units = "px") 

p1 = myqq(p_matrix[, "P0.001"]) + labs(tag = "A") + labs(title = "P+T(0.001)")
p2 = myqq(p_matrix[, "P0.05"]) + labs(tag = "B") + labs(title = "P+T(0.05)")
p3 = myqq(p_matrix[, "lassosum"]) + labs(tag = "C") + labs(title = "lassosum")
p4 = myqq(p_matrix[, "SDPR"]) + labs(tag = "D") + labs(title = "SDPR")
p5 = myqq(p_matrix[, "PRScs"]) + labs(tag = "E") + labs(title = "PRS-CS")
p6 = myqq(p_matrix[, "ACAT.res"]) + labs(tag = "F") + labs(title = "OTTERS")

grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

dev.off()

# Table 

sum_error = function(p_cut, mat){
  
  sum_vec = function(vec) return(sum(vec < p_cut, na.rm = T))
  
  return(apply(mat, 2, sum_vec))
  
}

Manu_dir = "~/YangFSSdata2/qdai/BS_TWAS/Scripts/Revision"
sapply(c(1e-2, 1e-4, 2.5 * 10^-6), sum_error, mat = p_matrix)
res.tab = sapply(c(1e-2, 1e-4, 2.5 * 10^-6), sum_error, mat = p_matrix)
res.tab = data.frame(t(res.tab / (500 * 2000)))
rownames(res.tab) = c(1e-2, 1e-4, 2.5 * 10^-6)
openxlsx::write.xlsx(res.tab,
                     file = file.path(Manu_dir, 
                                      'null_simu_results.xlsx'),
                     overwrite = T,
                     rowNames = T)

