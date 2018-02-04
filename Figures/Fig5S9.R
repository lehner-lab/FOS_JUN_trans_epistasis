library(parallel)


load("../Process_data/000-data/006-epistasis_thermo_model/trans_restricted.Rdata")
tr = d
load("../Process_data/000-data/006-epistasis_thermo_model/cis_restricted.Rdata")
ci = d


# p-value thresholds
seq_p = 10^seq(-5,0,length.out=100)

# total discoveries at this p-value threshold
# for positive and negative epistasis separately
td_pos_tr = cumsum(table(cut(tr$p[tr$mepi > 0],c(0,seq_p))))
td_neg_tr = cumsum(table(cut(tr$p[tr$mepi < 0],c(0,seq_p))))
td_pos_ci = cumsum(table(cut(ci$p[ci$mepi > 0],c(0,seq_p))))
td_neg_ci = cumsum(table(cut(ci$p[ci$mepi < 0],c(0,seq_p))))

# fdr at all p-value threshold
f_tr = sapply(seq_p,function(x){
    max(tr$fdr[tr$p < x])
})
f_ci = sapply(seq_p,function(x){
    max(ci$fdr[ci$p < x])
})

# proportion of true discoveries
prop_pos_tr = td_pos_tr * (1-f_tr) / nrow(tr)
prop_neg_tr = td_neg_tr * (1-f_tr) / nrow(tr)
prop_pos_ci = td_pos_ci * (1-f_ci) / nrow(ci)
prop_neg_ci = td_neg_ci * (1-f_ci) / nrow(ci)



# load from the 1000 matched libs
m_tr = do.call("cbind",mclapply(1:1000,function(x){
    load(paste0("../Process_data/000-data/009-fdr_matched_libs/trans",x,".Rdata"))
    prop
}))
m_ci = do.call("cbind",mclapply(1:1000,function(x){
    load(paste0("../Process_data/000-data/009-fdr_matched_libs/cis",x,".Rdata"))
    prop
}))


m_pos_mean_tr = rowMeans(m_tr[,rep(c(T,F),1000)])
m_neg_mean_tr = rowMeans(m_tr[,rep(c(F,T),1000)])
m_pos_err_tr = apply(m_tr[,rep(c(T,F),1000)],1,sd) / sqrt(1000) * 1.96
m_neg_err_tr = apply(m_tr[,rep(c(F,T),1000)],1,sd) / sqrt(1000) * 1.96
m_pos_mean_ci = rowMeans(m_ci[,rep(c(T,F),1000)])
m_neg_mean_ci = rowMeans(m_ci[,rep(c(F,T),1000)])
m_pos_err_ci = apply(m_ci[,rep(c(T,F),1000)],1,sd) / sqrt(1000) * 1.96
m_neg_err_ci = apply(m_ci[,rep(c(F,T),1000)],1,sd) / sqrt(1000) * 1.96


pdf("Fig5S9A.pdf")
plot(seq_p,m_pos_mean_tr,ylim=range(m_pos_mean_tr,m_pos_mean_ci,prop_pos_tr,prop_pos_ci),xlab="-log10(p-value)",ylab="Proportion of true positive genetic interactions")
segments(seq_p,m_pos_mean_tr-m_pos_err_tr,seq_p,m_pos_mean_tr+m_pos_err_tr)
points(seq_p,m_pos_mean_ci,col=2)
segments(seq_p,m_pos_mean_ci-m_pos_err_ci,seq_p,m_pos_mean_ci+m_pos_err_ci,col=2)
lines(seq_p,prop_pos_tr,lwd=2)
lines(seq_p,prop_pos_ci,col=2,lwd=2)
legend("topright",pch=c(1,1,NA,NA),lty=c(NA,NA,1,1),col=c(1,2,1,2),legend=c("matched trans","matched cis","original trans","original cis"))
dev.off()


pdf("Fig5S9B.pdf")
plot(seq_p,m_neg_mean_tr,ylim=range(m_neg_mean_tr,m_neg_mean_ci,prop_neg_tr,prop_neg_ci),xlab="-log10(p-value)",ylab="Proportion of true negative genetic interactions")
segments(seq_p,m_neg_mean_tr-m_neg_err_tr,seq_p,m_neg_mean_tr+m_neg_err_tr)
points(seq_p,m_neg_mean_ci,col=2)
segments(seq_p,m_neg_mean_ci-m_neg_err_ci,seq_p,m_neg_mean_ci+m_neg_err_ci,col=2)
lines(seq_p,prop_neg_tr,lwd=2)
lines(seq_p,prop_neg_ci,col=2,lwd=2)
legend("topright",pch=c(1,1,NA,NA),lty=c(NA,NA,1,1),col=c(1,2,1,2),legend=c("matched trans","matched cis","original trans","original cis"))
dev.off()


