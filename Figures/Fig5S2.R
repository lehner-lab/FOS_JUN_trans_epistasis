source("../000-functions.R")
library(gplots)

#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################

load("../Process_data/000-data/006-epistasis_thermo_model/trans_restricted.Rdata")
tr = d
load("../Process_data/000-data/006-epistasis_thermo_model/cis_restricted.Rdata")
ci = d


# define bins
nbin = 20
r = range(c(tr$s1_mppi,tr$s2_mppi,ci$s1_mppi,ci$s2_mppi))
s = seq(r[1],r[2],length.out=nbin+1)

# bin PPI scores
counts = cbind(
    table(cut(tr$s1_mppi,s, include.lowest = T)),
    table(cut(tr$s2_mppi,s, include.lowest = T)),
    table(cut(ci$s1_mppi,s, include.lowest = T)),
    table(cut(ci$s2_mppi,s, include.lowest = T))
)

# frequencies
ori = apply(counts,2,function(x){x/sum(x)})

freq = mclapply(1:1000,function(x){
    load(paste0("../Process_data/000-data/007-matched_libs/trans",x,".Rdata"))
    load(paste0("../Process_data/000-data/007-matched_libs/cis",x,".Rdata"))
    cbind(
        table(cut(tr$s1_mppi,s, include.lowest = T)),
        table(cut(tr$s2_mppi,s, include.lowest = T)),
        table(cut(ci$s1_mppi,s, include.lowest = T)),
        table(cut(ci$s2_mppi,s, include.lowest = T))
    ) / nrow(tr) # since they're matched, tr and ci have the same size
},mc.cores=8)
freq = simplify2array(freq)


# mean and confidnec einterval over the 1000 matched-libraries
fr_mean = sapply(1:4,function(x){
    rowMeans(freq[,x,])
})
fr_conf = sapply(1:4,function(x){
    apply(freq[,x,],1,sd) / sqrt(1000) * 1.96
})



pdf("Fig5S2A.pdf")
yl = range(fr_mean+fr_conf,ori)
ticks = cbind(round(s,digits=2),seq(0,40,2))

b = barplot2(t(fr_mean[,c(1,3)]), beside=T, space=c(0,0), plot.ci=T, ci.l=t(fr_mean[,c(1,3)]) - t(fr_conf[,c(1,3)]), ci.u=t(fr_mean[,c(1,3)]) + t(fr_conf[,c(1,3)]),legend.text=c("trans","cis"), names.arg=NULL, ylim=yl, ylab="Density", xlab="PPI score mutant 1")
axis(1,at=ticks[seq(1,21,5),2],labels=ticks[seq(1,21,5),1])
lines(colMeans(b), ori[,1], lwd=4, col=2)
lines(colMeans(b), ori[,3], lwd=4, col=7)

b = barplot2(t(fr_mean[,c(2,4)]), beside=T, space=c(0,0), plot.ci=T, ci.l=t(fr_mean[,c(2,4)]) - t(fr_conf[,c(2,4)]), ci.u=t(fr_mean[,c(2,4)]) + t(fr_conf[,c(2,4)]),legend.text=c("trans","cis"), names.arg=NULL, ylim=yl, ylab="Density", xlab="PPI score mutant 2")
axis(1,at=ticks[seq(1,21,5),2],labels=ticks[seq(1,21,5),1])
lines(colMeans(b), ori[,2], lwd=4, col=2)
lines(colMeans(b), ori[,4], lwd=4, col=7)
dev.off()



#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################

r = range(c(tr$do_mppi,ci$do_mppi))
s = seq(1, r[2], 0.01)

freq = mclapply(1:1000,function(x){
    load(paste0("../Process_data/000-data/007-matched_libs/trans",x,".Rdata"))
    load(paste0("../Process_data/000-data/007-matched_libs/cis",x,".Rdata"))
    cbind(
        rev(cumsum(rev(table(cut(tr$do_mppi,s, include.lowest = T))))),
        rev(cumsum(rev(table(cut(ci$do_mppi,s, include.lowest = T)))))
    ) / nrow(tr) # since they're matched, tr and ci have the same size
},mc.cores=8)
freq = simplify2array(freq)


fr_mean = sapply(1:2,function(x){
    rowMeans(freq[,x,])
})
fr_conf = sapply(1:2,function(x){
    apply(freq[,x,],1,sd) / sqrt(1000) * 1.96
})

pdf("Fig5S2B.pdf")
# confidence intervals are not added they are in the order of magnitude of 1e-5
plot(s[1:(length(s)-1)], fr_mean[,1], type="l", lwd=3, ylim=range(fr_mean), xlab="threshold above which double mutants are classified as stronger", ylab="proportion")
lines(s[1:(length(s)-1)], fr_mean[,2], lwd=3, col=2)
legend("topright",lwd=3,col=c(1,2),legend=c("trans","cis"))
dev.off()





#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################

s = seq(r[1], 0.8, 0.01)

freq = mclapply(1:1000,function(x){
    load(paste0("../Process_data/000-data/007-matched_libs/trans",x,".Rdata"))
    load(paste0("../Process_data/000-data/007-matched_libs/cis",x,".Rdata"))
    cbind(
        cumsum(table(cut(tr$do_mppi,s, include.lowest = T))),
        cumsum(table(cut(ci$do_mppi,s, include.lowest = T)))
    ) / nrow(tr) # since they're matched, tr and ci have the same size
},mc.cores=8)
freq = simplify2array(freq)


fr_mean = sapply(1:2,function(x){
    rowMeans(freq[,x,])
})
fr_conf = sapply(1:2,function(x){
    apply(freq[,x,],1,sd) / sqrt(1000) * 1.96
})

pdf("Fig5S2C.pdf")
# confidence intervals are not added they are in the order of magnitude of 1e-5
plot(s[1:(length(s)-1)], fr_mean[,1], type="l", lwd=3, ylim=range(fr_mean), xlab="threshold under which double mutants are classified as severly detrimental", ylab="proportion")
lines(s[1:(length(s)-1)], fr_mean[,2], lwd=3, col=2)
legend("topleft",lwd=3,col=c(1,2),legend=c("trans","cis"))
dev.off()



#####################################################################################################################
### Panel D  ########################################################################################################
#####################################################################################################################

s = seq(0.5, 0.8, 0.01)

freq = mclapply(1:1000,function(x){
    load(paste0("../Process_data/000-data/007-matched_libs/trans",x,".Rdata"))
    load(paste0("../Process_data/000-data/007-matched_libs/cis",x,".Rdata"))
    cbind(
        rev(cumsum(rev(table(cut(tr$do_mppi[tr$do_mppi < 0.92],s, include.lowest = T))))),
        rev(cumsum(rev(table(cut(ci$do_mppi[ci$do_mppi < 0.92],s, include.lowest = T)))))
    ) / nrow(tr) # since they're matched, tr and ci have the same size
},mc.cores=8)
freq = simplify2array(freq)


fr_mean = sapply(1:2,function(x){
    rowMeans(freq[,x,])
})
fr_conf = sapply(1:2,function(x){
    apply(freq[,x,],1,sd) / sqrt(1000) * 1.96
})

pdf("Fig5S2D.pdf")
# confidence intervals are not added they are in the order of magnitude of 1e-5
plot(s[1:(length(s)-1)], fr_mean[,1], type="l", lwd=3, ylim=range(fr_mean), xlab="threshold above which double mutants are classified as intermediate (top threshold = 0.92)", ylab="proportion")
lines(s[1:(length(s)-1)], fr_mean[,2], lwd=3, col=2)
legend("topright",lwd=3,col=c(1,2),legend=c("trans","cis"))
dev.off()



#####################################################################################################################
### Panel E  ########################################################################################################
#####################################################################################################################

s = seq(0.64, 1, 0.01)

freq = mclapply(1:1000,function(x){
    load(paste0("../Process_data/000-data/007-matched_libs/trans",x,".Rdata"))
    load(paste0("../Process_data/000-data/007-matched_libs/cis",x,".Rdata"))
    cbind(
        cumsum(table(cut(tr$do_mppi[tr$do_mppi > 0.64],s, include.lowest = T))),
        cumsum(table(cut(ci$do_mppi[ci$do_mppi > 0.64],s, include.lowest = T)))
    ) / nrow(tr) # since they're matched, tr and ci have the same size
},mc.cores=8)
freq = simplify2array(freq)


fr_mean = sapply(1:2,function(x){
    rowMeans(freq[,x,])
})
fr_conf = sapply(1:2,function(x){
    apply(freq[,x,],1,sd) / sqrt(1000) * 1.96
})

pdf("Fig5S2E.pdf")
# confidence intervals are not added they are in the order of magnitude of 1e-5
plot(s[1:(length(s)-1)], fr_mean[,1], type="l", lwd=3, ylim=range(fr_mean), xlab="threshold under which double mutants are classified as intermediate (bottom threshold = 0.64)", ylab="proportion")
lines(s[1:(length(s)-1)], fr_mean[,2], lwd=3, col=2)
legend("topleft",lwd=3,col=c(1,2),legend=c("trans","cis"))
dev.off()



