library(ggplot2)
library(gridExtra)
source("../000-functions.R")

#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################


# load PPi scores calculated at different read quakity threshold
# this is data with read counts filtered since, as will be shown in panel B, low read counts give biased estimates of PPI scores
load("../Process_data/000-data/002-PPIscores/trans_Q20_readCountFiltered.Rdata")
q20 = d
load("../Process_data/000-data/002-PPIscores/trans_Q30_readCountFiltered.Rdata")
q30 = d
load("../Process_data/000-data/002-PPIscores/trans_Q35_readCountFiltered.Rdata")
q35 = d
load("../Process_data/000-data/002-PPIscores/trans_mutQualFiltered_Q20_readCountFiltered.Rdata")
q20f = d


# merge with the PPI scores that will be used for the rest of the analyses, i.e. filtered at average Phred score of 20
plot1 = merge(q20,q30,by=c("id1","id2"))
plot2 = merge(q20,q35,by=c("id1","id2"))
plot3 = merge(q20,q20f,by=c("id1","id2"))
cor1 = round(cor(plot1$mppi.x, plot1$mppi.y),digits=3)
cor2 = round(cor(plot2$mppi.x, plot2$mppi.y),digits=3)
cor3 = round(cor(plot3$mppi.x, plot3$mppi.y),digits=3)

# plot
pdf("Fig1S1A.pdf")
par(mfrow=c(2,2))
smoothScatter(plot1$mppi.x, plot1$mppi.y, xlab = "PPI scores", ylab = "PPI scores with read quality filtered at Q30")
text(0.2,1.1,paste0("R=",cor1),pos=4)
smoothScatter(plot2$mppi.x, plot2$mppi.y, xlab = "PPI scores", ylab = "PPI scores with read quality filtered at Q35")
text(0.2,1.1,paste0("R=",cor2),pos=4)
smoothScatter(plot3$mppi.x, plot3$mppi.y, xlab = "PPI scores", ylab = "PPI scores with low quality mutations filtered out")
text(0.2,1.1,paste0("R=",cor3),pos=4)
dev.off()




#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################

load("../Process_data/000-data/002-PPIscores/trans_Q20.Rdata")

g1= ggplot(data = d, aes(i1, ppi1)) +
  stat_density2d(aes(fill = ..density..^0.25), geom = "tile", contour = FALSE, n = 200) +
  scale_x_continuous(trans = 'log10') + 
  scale_fill_continuous(low = "white", high = "#08306B") +
  xlab("log10(read count) replicate 1") + 
  ylab("PPI scores replicate 1") + 
  theme_classic() + 
  geom_vline(xintercept=10)

g2 = ggplot(data = d, aes(ppi1)) +
  geom_histogram(bins = 40) +
  coord_flip() + 
  xlab("PPI scores replicate 1") + 
  scale_y_reverse()


pdf("Fig1S1B.pdf")
grid.arrange(g2, g1, ncol=2, nrow=1, widths=c(1, 4))
dev.off()



#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################


# correlation between replicates, using the q20 library and filtered based on input and output read counts
pdf("Fig1S1C.pdf")
pairs(q20[,paste0("ppi",1:3)],pch=".", upper.panel=panel.cor.pearson, diag.panel=panel.hist40, lower.panel=panel.smooth)
dev.off()