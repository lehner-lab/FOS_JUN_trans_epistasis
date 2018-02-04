library(RColorBrewer)

#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################


load("../Process_data/000-data/002-PPIscores/trans_Q20_readCountFiltered.Rdata")

# split into doubles and singles
do = d[d$pos1 != 0 & d$pos2 != 0,]
si = d[d$pos1 == 0 | d$pos2 == 0,]

cor_do = round(cor(do$ppi1, do$ppi2), digits=3)
cor_si = round(cor(si$ppi1, si$ppi2), digits=3)

pdf("Fig1C.pdf")
smoothScatter(do$ppi1, do$ppi2, xlab = "PPI scores replicate 1", ylab = "PPI scores replicate 2")
points(si$ppi1, si$ppi2, pch=16, col=brewer.pal(9,"Set3")[6], cex=0.6)
abline(0,1)
legend("topleft", legend = c(paste("singles, R =", cor_si), paste("doubles, R =", cor_do)), fill=c(NA,blues9[5]), pch=c(16,NA), col=c(brewer.pal(9,"Set3")[6],NA), border=NA)
dev.off()





#####################################################################################################################
### Panel D  ########################################################################################################
#####################################################################################################################

# load ODs
od = read.delim("000-additional_data/001-small_scale_conf_ODs.txt")

# remove first reading because it is generally off for some reasons
od = od[,names(od) != "X0",]

# keep only MTX condition
od = od[od$cond == "MTX",]

# subtract blank and convert ODs to log
empty = od[od$prot == "empty",]
od = od[od$prot != "empty",]
od[,7:ncol(od)] = log(od[,7:ncol(od)] - colMeans(empty[,7:ncol(empty)]))

# time vector (in hours)
t = seq(0.25,0.25*(ncol(od)-6),0.25)

# identify the time point at which the WT has its max growth rate
wt = as.matrix(od[od$prot == "wt",7:ncol(od)])
# time sliding window length for the regression. Has to be an odd number
w.l = 31
# for each replicate
wt.gr = sapply(1:2,function(x){
  # for each sliding window
  do.call("c",mclapply(1:(length(t)-w.l+1),function(start){
    # do a linear regression to get the slope, i.e. the growth rate
    lm(wt[x,start:(start+w.l-1)] ~ t[start:(start+w.l-1)])$coefficients[2]
  },mc.cores=8))
})
max.wt.gr = apply(wt.gr,2,max)
idx.max.gr = mean(apply(wt.gr,2,function(x){
  which(x == max(x)) + floor(w.l/2)
}))



# compute the equivalent of the PPI score for each variant
# this is done by computing the number of generations since the beginning of the experiment until the time at which the wt reached it's max growth rate
# the ODs are averaged around the time points to average out measurement errors
window.av = 5
# time points around max wt grwoth
output = (idx.max.gr-floor(window.av/2)):(idx.max.gr+floor(window.av/2))

od2 = split(od[,7:ncol(od)],list(od$prot,od$mut), drop=T)
# for each mutants
gen = mclapply(od2,function(x){
  # for each replicate
  do.call("c",lapply(1:nrow(x),function(i){
    x = as.numeric(as.character(x[i,]))
    # number of generation grown = log(ratio of ODs). Data are already logged so it is the difference of log(ODs)
    (mean(x[output]) - mean(x[1:window.av]))
  }))
},mc.cores=8)

wt.gen = gen[["wt.0gwt"]]

ppi = do.call("rbind",lapply(gen,function(x){
  # average ppi score
  m = mean(x) / mean(wt.gen)
  # 95% confidence interval
  s = sd(x / mean(wt.gen)) / sqrt(length(x)) * 1.96
  c(m,s)
}))
ppi = as.data.frame(ppi)
names(ppi) = c("meanPPI","err")


# compare to deepPCA PPI scores
ids = do.call("rbind",strsplit(row.names(ppi),"\\."))
ppi$id1[ids[,1] == "FOS"] = ids[ids[,1] == "FOS",2]
ppi$id2[ids[,1] == "JUN"] = ids[ids[,1] == "JUN",2]
ppi$id1[ids[,1] != "FOS"] = "0gwt"
ppi$id2[ids[,1] != "JUN"] = "0gwt"


d = merge(ppi,si,by=c("id1","id2"),all.x=T)
# confidence interval for the deepPCA PPI scores
d$errLS = apply(d[,paste0("ppi",1:3)],1,sd) / sqrt(3) * 1.96

# WT had been filtered out, so add its PPI score, which is 1 by definition
d$mppi[is.na(d$mppi)] = 1
d$errLS[is.na(d$errLS)] = 0

pdf("Fig1D.pdf")
co = round(cor(d$mppi,d$meanPPI),digits=3)
plot(d$mppi,d$meanPPI,pch=16,cex=1.3,col="darkblue",ylab="growth rate relative to wt (TECAN measurements)",xlab="PPI score",main = paste("R =",co), xlim=c(min(d$mppi-d$errLS),max(d$mppi+d$errLS)), ylim=c(min(d$meanPPI-d$err),max(d$meanPPI+d$err)))
l = lm(d$meanPPI ~ d$mppi)
abline(l,col="grey80",lty=2)
segments(d$mppi-d$errLS,d$meanPPI,d$mppi+d$errLS,d$meanPPI,col="darkblue",lwd=2)
segments(d$mppi,d$meanPPI-d$err,d$mppi,d$meanPPI+d$err,col="darkblue",lwd=2)
dev.off()
