source("../000-functions.R")

# load multiplicative epistasis data
load("../Process_data/000-data/003-mult_epi/trans.Rdata")

pdf("Fig3S1.pdf")
pairs(d[,paste0("epi",1:3)], upper.panel = panel.cor.pearson, lower.panel = panel.smooth, diag.panel = panel.hist40)
dev.off()