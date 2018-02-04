source("../000-functions.R")


#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################

load("../Process_data/000-data/006-epistasis_thermo_model/trans_restricted.Rdata")
tr = d
load("../Process_data/000-data/006-epistasis_thermo_model/cis_restricted.Rdata")
ci = d


pdf("Fig5S4A.pdf")
pairs(tr[,paste0("epi",1:3)],pch=".", upper.panel=panel.cor.pearson, diag.panel=panel.hist40, lower.panel=panel.smooth, main="trans")
pairs(ci[,paste0("epi",1:3)],pch=".", upper.panel=panel.cor.pearson, diag.panel=panel.hist40, lower.panel=panel.smooth, main="cis")
dev.off()


#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################



pdf("Fig5S4B.pdf")

# trans
p_thresh = min(tr$p[!tr$hit])
plot(tr$mepi, -log10(tr$p), xlab="Genetic interaction score", ylab="-log10(p-value)", pch=".", main="trans")
abline(h=-log10(p_thresh), lwd=2, col=2)
abline(v=c(-0.1,0.1), col=2, lwd=2)
text(-0.45,4.5,"negative genetic interactions", pos=4)
text(0.45,4.5,"positive genetic interactions", pos=2)
text(-0.48, -log10(p_thresh)-0.4, paste0("p < ", round(p_thresh,digits=6),"\nFDR < 20%"), col=2, pos=4)

# cis
p_thresh = min(ci$p[!ci$hit])
plot(ci$mepi, -log10(ci$p), xlab="Genetic interaction score", ylab="-log10(p-value)", pch=".", main="cis")
abline(h=-log10(p_thresh), lwd=2, col=2)
abline(v=c(-0.1,0.1), col=2, lwd=2)
text(-0.6,4.5,"negative genetic interactions", pos=4)
text(0.4,4.5,"positive genetic interactions", pos=2)
text(-0.6, -log10(p_thresh)-0.4, paste0("p < ", round(p_thresh,digits=6),"\nFDR < 20%"), col=2, pos=4)
dev.off()

