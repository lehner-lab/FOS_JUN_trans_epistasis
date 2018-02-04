library(parallel)
library(aqfig)


load("../Process_data/000-data/006-epistasis_thermo_model/trans.Rdata")




# epistasis thresholds
s = seq(0,0.6,0.05)

unique.pos = sort(unique(substr(d$id1,1,2)))

r = do.call("rbind",lapply(as.list(s),function(epi.thresh){

    # hits at this threshold
    hp = factor(d$mepi > epi.thresh & d$hit, levels=c(T,F))
    hn = factor(d$mepi < -epi.thresh & d$hit, levels=c(T,F))
    
    # for each pair of position
    r = do.call("rbind",lapply(as.list(1:32),function(x){
        do.call("rbind",mclapply(as.list(1:32),function(y){
            # this pos pair
            pos = factor(d$pos1 == x & d$pos2 == y, levels=c(T,F))
            # fisher test
            fp = fisher.test(table(pos,hp))
            fn = fisher.test(table(pos,hn))
            c(fp$estimate,fp$p.value,fn$estimate,fn$p.value)
        },mc.cores=8))
    }))
    
    # FDR is calculated at each threshold independently since it to compare the results if we had done the analyses at this threshold
    r = cbind(r,p.adjust(r[,2],method="fdr"))
    r = cbind(r,p.adjust(r[,4],method="fdr"))
    r = as.data.frame(r)
    r = cbind(epi.thresh, rep(unique.pos,each=32),rep(unique.pos,32),r)
    
}))



load("Fig4_data_forS2.Rdata")

enriched_pairs_pos = paste(r2[r2[,6] == 2,1],r2[r2[,6] == 2,2])
enriched_pairs_neg = paste(r2[r2[,6] == 1,1],r2[r2[,6] == 1,2])

# keep only enriched pairs
r.pos = r[paste(r[,2],r[,3]) %in% enriched_pairs_pos,]
r2.pos = matrix(log2(r.pos[,4]), ncol=length(s))
colnames(r2.pos) = s
rownames(r2.pos) = paste(r.pos[r.pos[,1] == 0,2],r.pos[r.pos[,1] == 0,3])
r2.pos = t(r2.pos)

r.neg = r[paste(r[,2],r[,3]) %in% enriched_pairs_neg,]
r2.neg = matrix(log2(r.neg[,6]), ncol=length(s))
colnames(r2.neg) = s
rownames(r2.neg) = paste(r.neg[r.neg[,1] == 0,2],r.neg[r.neg[,1] == 0,3])
r2.neg = t(r2.neg)
# add fake columns just to keep the size of the graphics the same
r2.neg = cbind(r2.neg,matrix(NA,nrow=nrow(r2.neg),ncol=ncol(r2.pos)-ncol(r2.neg)))

# compute FDR and identify non-significant enrichments
r2.pos.fdr = t(matrix(r.pos[,5],nrow=nrow(r2[r2[,6] == 2,])))
r2.neg.fdr = t(matrix(r.neg[,7],nrow=nrow(r2[r2[,6] == 1,])))
r2.neg.fdr = cbind(r2.neg.fdr,matrix(NA,nrow=nrow(r2.neg.fdr),ncol=ncol(r2.pos.fdr)-ncol(r2.neg.fdr)))
r2.pos.fdr = which(r2.pos.fdr >= 0.1 & is.finite(r2.pos), arr.ind=T)
r2.neg.fdr = which(r2.neg.fdr >= 0.1 & is.finite(r2.neg), arr.ind=T)


pdf("Fig4S2.pdf",height=14)
par(mfrow=c(2,1),mar=c(5,5,5,5))
m = ceiling(max(c(r2.pos[is.finite(r2.pos)],r2.neg[is.finite(r2.neg)])))
zl = c(0,m)
image(1:nrow(r2.pos),1:ncol(r2.pos),r2.pos,axes=F,xlab="epistasis magnitude threshold",ylab="enriched pairs",col=heat.colors(m-1),zlim=zl)
axis(1,at=1:nrow(r2.pos),labels=rownames(r2.pos))
axis(2,at=1:ncol(r2.pos),labels=colnames(r2.pos),las=2)
points(r2.pos.fdr[,1],r2.pos.fdr[,2],pch=16)

image(1:nrow(r2.neg),1:ncol(r2.neg),r2.neg,axes=F,xlab="epistasis magnitude threshold",ylab="enriched pairs",col=heat.colors(m),zlim=zl)
axis(1,at=1:nrow(r2.neg),labels=rownames(r2.neg))
axis(2,at=1:ncol(r2.neg),labels=colnames(r2.neg),las=2)
points(r2.neg.fdr[,1],r2.neg.fdr[,2],pch=16)

vertical.image.legend(zlim=zl,col=heat.colors(m))
dev.off()