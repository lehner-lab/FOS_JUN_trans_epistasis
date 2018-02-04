source("../000-functions.R")


#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################

load("../Process_data/000-data/002-PPIscores/cis_Q20_restricted.Rdata")
ci = d

pdf("Fig5S1A.pdf")
pairs(ci[,paste0("ppi",1:3)],pch=".", upper.panel=panel.cor.pearson, diag.panel=panel.hist40, lower.panel=panel.smooth)
dev.off()


#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################

load("../Process_data/000-data/002-PPIscores/trans_Q20_restricted.Rdata")
tr = d

# keep only singles
tr = tr[tr$pos2 == 0,]
ci = ci[ci$pos1 == 0,]

d = merge(tr,ci,by.x="id1",by.y="id2")

co = round(cor(d$mppi.x, d$mppi.y),digits=2)

pdf("Fig5S1B.pdf")
plot(d$mppi.x, d$mppi.y, xlab="PPI scores trans", ylab="PPI scores cis")
abline(0,1)
text(0.5,1,paste0("R=",co))
dev.off()
