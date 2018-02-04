library(parallel)
library(RColorBrewer)
library(aqfig)

load("../Process_data/000-data/006-epistasis_thermo_model/trans_restricted.Rdata")
tr = d
load("../Process_data/000-data/006-epistasis_thermo_model/cis_restricted.Rdata")
ci = d


#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################


# trans first
d = tr

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



load("Fig5S6_data_forS8_trans.Rdata")

enriched_pairs_pos = paste(r2.trans[r2.trans[,6] == 2,1],r2.trans[r2.trans[,6] == 2,2])
enriched_pairs_neg = paste(r2.trans[r2.trans[,6] == 1,1],r2.trans[r2.trans[,6] == 1,2])

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
r2.pos.fdr = t(matrix(r.pos[,5],nrow=nrow(r2.trans[r2.trans[,6] == 2,])))
r2.neg.fdr = t(matrix(r.neg[,7],nrow=nrow(r2.trans[r2.trans[,6] == 1,])))
r2.neg.fdr = cbind(r2.neg.fdr,matrix(NA,nrow=nrow(r2.neg.fdr),ncol=ncol(r2.pos.fdr)-ncol(r2.neg.fdr)))
r2.pos.fdr = which(r2.pos.fdr >= 0.1 & is.finite(r2.pos), arr.ind=T)
r2.neg.fdr = which(r2.neg.fdr >= 0.1 & is.finite(r2.neg), arr.ind=T)


pdf("Fig5S8A_trans.pdf",height=14)
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






# cis
d = ci

# epistasis thresholds
s = seq(0,0.6,0.05)

unique.pos = sort(unique(substr(c(d$id1,d$id2),1,2)))
names(unique.pos) = 1:32

r = do.call("rbind",lapply(as.list(s),function(epi.thresh){
    
    # hits at this threshold
    hp = factor(d$mepi > epi.thresh & d$hit, levels=c(T,F))
    hn = factor(d$mepi < -epi.thresh & d$hit, levels=c(T,F))
    
    # for each pair of position
    r = do.call("rbind",lapply(as.list(1:31),function(x){
        do.call("rbind",mclapply(as.list((x+1):32),function(y){
            # this pos pair
            pos = factor(d$pos1 == x & d$pos2 == y, levels=c(T,F))
            # fisher test
            fp = fisher.test(table(pos,hp))
            fn = fisher.test(table(pos,hn))
            c(x,y,fp$estimate,fp$p.value,fn$estimate,fn$p.value)
        },mc.cores=8))
    }))
    
    # FDR is calculated at each threshold independently since it to compare the results if we had done the analyses at this threshold
    r = cbind(r,p.adjust(r[,4],method="fdr"))
    r = cbind(r,p.adjust(r[,6],method="fdr"))
    r = as.data.frame(r)
    r = cbind(epi.thresh, r)
    
}))

r[,2] = unique.pos[r[,2]]
r[,3] = unique.pos[r[,3]]



load("Fig5S6_data_forS8_cis.Rdata")

enriched_pairs_pos = paste(r2.cis[r2.cis[,6] == 2,1],r2.cis[r2.cis[,6] == 2,2])
enriched_pairs_neg = paste(r2.cis[r2.cis[,6] == 1,1],r2.cis[r2.cis[,6] == 1,2])

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
r2.pos.fdr = t(matrix(r.pos[,5],nrow=nrow(r2.cis[r2.cis[,6] == 2,])))
r2.neg.fdr = t(matrix(r.neg[,7],nrow=nrow(r2.cis[r2.cis[,6] == 1,])))
r2.neg.fdr = cbind(r2.neg.fdr,matrix(NA,nrow=nrow(r2.neg.fdr),ncol=ncol(r2.pos.fdr)-ncol(r2.neg.fdr)))
r2.pos.fdr = which(r2.pos.fdr >= 0.1 & is.finite(r2.pos), arr.ind=T)
r2.neg.fdr = which(r2.neg.fdr >= 0.1 & is.finite(r2.neg), arr.ind=T)


pdf("Fig5S8A_cis.pdf",height=14)
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





#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################


# trans first
d = tr


# type of position
d$type1[substr(d$id1,2,2) %in% c("a","d")] = "core"
d$type1[substr(d$id1,2,2) %in% c("e","g")] = "salt b."
d$type1[substr(d$id1,2,2) %in% c("b","c","f")] = "far"
d$type2[substr(d$id2,2,2) %in% c("a","d")] = "core"
d$type2[substr(d$id2,2,2) %in% c("e","g")] = "salt b."
d$type2[substr(d$id2,2,2) %in% c("b","c","f")] = "far"


d$type1 = factor(d$type1)
d$type2 = factor(d$type2)

# epistasis threshold
s = seq(0,0.5,0.01)

# for negative epistasis
# six first columns are log odds ratios, 6 last are FDR
r.n = do.call("rbind",mclapply(as.list(s),function(epi.thresh){
    
    # number of hits for each pair of types
    t1 = table(d[d$mepi < -epi.thresh & d$hit,c("type1","type2")])
    t2 = table(d[!(d$mepi < -epi.thresh & d$hit),c("type1","type2")])
    
    # build contingency tables and do Fisher test
    r = do.call("rbind",lapply(as.list(1:3),function(i){
        do.call("rbind",lapply(as.list(i:3),function(j){
            if (i == j){
                m = matrix(c(t1[i,j],sum(t1)-t1[i,j],t2[i,j],sum(t2)-t2[i,j]),nrow=2) # rows = type / pas type, cols = hit / pas hit
            }else{
                m = matrix(c(t1[i,j]+t1[j,i],sum(t1)-t1[i,j]-t1[j,i],t2[i,j]+t2[j,i],sum(t2)-t2[i,j]-t2[j,i]),nrow=2)  
            }
            f = fisher.test(m)
            c(log(f$estimate),f$p.value)
        }))
    }))
    
    # FDR
    r[,2] = p.adjust(r[,2],method="fdr")
    
    as.vector(r)
}))
r.n[r.n == -Inf] = min(r.n[,1:6][is.finite(r.n[,1:6])])
r.n[r.n == Inf] = max(r.n[,1:6][is.finite(r.n[,1:6])])


# same for positive epistasis
r.p = do.call("rbind",mclapply(as.list(s),function(epi.thresh){
    t1 = table(d[d$mepi > epi.thresh & d$hit,c("type1","type2")])
    t2 = table(d[!(d$mepi > epi.thresh & d$hit),c("type1","type2")])
    
    r = do.call("rbind",lapply(as.list(1:3),function(i){
        do.call("rbind",lapply(as.list(i:3),function(j){
            if (i == j){
                m = matrix(c(t1[i,j],sum(t1)-t1[i,j],t2[i,j],sum(t2)-t2[i,j]),nrow=2) # rows = type / pas type, cols = hit / pas hit
            }else{
                m = matrix(c(t1[i,j]+t1[j,i],sum(t1)-t1[i,j]-t1[j,i],t2[i,j]+t2[j,i],sum(t2)-t2[i,j]-t2[j,i]),nrow=2)  
            }
            f = fisher.test(m)
            c(log(f$estimate),f$p.value)
        }))
    }))
    
    # FDR
    r[,2] = p.adjust(r[,2],method="fdr")
    
    as.vector(r)
}))
r.p[r.p == -Inf] = min(r.p[,1:6][is.finite(r.p[,1:6])])
r.p[r.p == Inf] = max(r.p[,1:6][is.finite(r.p[,1:6])])





pdf("Fig5S8B_trans.pdf")

# positive epistasis
# plot the enrichment lines for each pair of position type
plot(1,type="n",xlim=range(s),ylim=range(r.p[,1:6]),xlab="genetic interaction magnitude threshold",ylab="log10(odds ratio)", main="positive genetic interactions")
lapply(as.list(1:6),function(x){
    lines(s,r.p[,x],col=x,lwd=2)
})
legend("bottomleft",col=1:6,lwd=2,legend=c("core x core","far x core","salt x core","far x far","far x salt","salt x salt"))

# identify significant enrichments
# first remove FDR > 10%
for(i in 1:6){
    for(j in 1:length(s)){
        if(r.p[j,i+6] > 0.1){
            r.p[j,i] = NA
        }
    }
}
# then add a cross for significant ones
lapply(as.list(1:6),function(x){
    points(s,r.p[,x],col=x,pch=4)
})


# negative epistasis
# plot the enrichment lines for each pair of position type
plot(1,type="n",xlim=range(s),ylim=range(r.n[,1:6]),xlab="genetic interaction magnitude threshold",ylab="log10(odds ratio)", main="negative genetic interactions")
lapply(as.list(1:6),function(x){
    lines(s,r.n[,x],col=x,lwd=2)
})
legend("bottomleft",col=1:6,lwd=2,legend=c("core x core","far x core","salt x core","far x far","far x salt","salt x salt"))

# identify significant enrichments
# first remove FDR > 10%
for(i in 1:6){
    for(j in 1:length(s)){
        if(r.n[j,i+6] > 0.1){
            r.n[j,i] = NA
        }
    }
}
# then add a cross for significant ones
lapply(as.list(1:6),function(x){
    points(s,r.n[,x],col=x,pch=4)
})
dev.off()






# cis
d = ci


# type of position
d$type1[substr(d$id1,2,2) %in% c("a","d")] = "core"
d$type1[substr(d$id1,2,2) %in% c("e","g")] = "salt b."
d$type1[substr(d$id1,2,2) %in% c("b","c","f")] = "far"
d$type2[substr(d$id2,2,2) %in% c("a","d")] = "core"
d$type2[substr(d$id2,2,2) %in% c("e","g")] = "salt b."
d$type2[substr(d$id2,2,2) %in% c("b","c","f")] = "far"


d$type1 = factor(d$type1)
d$type2 = factor(d$type2)

# epistasis threshold
s = seq(0,0.5,0.01)

# for negative epistasis
# six first columns are log odds ratios, 6 last are FDR
r.n = do.call("rbind",mclapply(as.list(s),function(epi.thresh){
    
    # number of hits for each pair of types
    t1 = table(d[d$mepi < -epi.thresh & d$hit,c("type1","type2")])
    t2 = table(d[!(d$mepi < -epi.thresh & d$hit),c("type1","type2")])
    
    # build contingency tables and do Fisher test
    r = do.call("rbind",lapply(as.list(1:3),function(i){
        do.call("rbind",lapply(as.list(i:3),function(j){
            if (i == j){
                m = matrix(c(t1[i,j],sum(t1)-t1[i,j],t2[i,j],sum(t2)-t2[i,j]),nrow=2) # rows = type / pas type, cols = hit / pas hit
            }else{
                m = matrix(c(t1[i,j]+t1[j,i],sum(t1)-t1[i,j]-t1[j,i],t2[i,j]+t2[j,i],sum(t2)-t2[i,j]-t2[j,i]),nrow=2)  
            }
            f = fisher.test(m)
            c(log(f$estimate),f$p.value)
        }))
    }))
    
    # FDR
    r[,2] = p.adjust(r[,2],method="fdr")
    
    as.vector(r)
}))
r.n[r.n == -Inf] = min(r.n[,1:6][is.finite(r.n[,1:6])])
r.n[r.n == Inf] = max(r.n[,1:6][is.finite(r.n[,1:6])])


# same for positive epistasis
r.p = do.call("rbind",mclapply(as.list(s),function(epi.thresh){
    t1 = table(d[d$mepi > epi.thresh & d$hit,c("type1","type2")])
    t2 = table(d[!(d$mepi > epi.thresh & d$hit),c("type1","type2")])
    
    r = do.call("rbind",lapply(as.list(1:3),function(i){
        do.call("rbind",lapply(as.list(i:3),function(j){
            if (i == j){
                m = matrix(c(t1[i,j],sum(t1)-t1[i,j],t2[i,j],sum(t2)-t2[i,j]),nrow=2) # rows = type / pas type, cols = hit / pas hit
            }else{
                m = matrix(c(t1[i,j]+t1[j,i],sum(t1)-t1[i,j]-t1[j,i],t2[i,j]+t2[j,i],sum(t2)-t2[i,j]-t2[j,i]),nrow=2)  
            }
            f = fisher.test(m)
            c(log(f$estimate),f$p.value)
        }))
    }))
    
    # FDR
    r[,2] = p.adjust(r[,2],method="fdr")
    
    as.vector(r)
}))
r.p[r.p == -Inf] = min(r.p[,1:6][is.finite(r.p[,1:6])])
r.p[r.p == Inf] = max(r.p[,1:6][is.finite(r.p[,1:6])])





pdf("Fig5S8B_cis.pdf")

# positive epistasis
# plot the enrichment lines for each pair of position type
plot(1,type="n",xlim=range(s),ylim=range(r.p[,1:6]),xlab="genetic interaction magnitude threshold",ylab="log10(odds ratio)", main="positive genetic interactions")
lapply(as.list(1:6),function(x){
    lines(s,r.p[,x],col=x,lwd=2)
})
legend("bottomleft",col=1:6,lwd=2,legend=c("core x core","far x core","salt x core","far x far","far x salt","salt x salt"))

# identify significant enrichments
# first remove FDR > 10%
for(i in 1:6){
    for(j in 1:length(s)){
        if(r.p[j,i+6] > 0.1){
            r.p[j,i] = NA
        }
    }
}
# then add a cross for significant ones
lapply(as.list(1:6),function(x){
    points(s,r.p[,x],col=x,pch=4)
})


# negative epistasis
# plot the enrichment lines for each pair of position type
plot(1,type="n",xlim=range(s),ylim=range(r.n[,1:6]),xlab="genetic interaction magnitude threshold",ylab="log10(odds ratio)", main="negative genetic interactions")
lapply(as.list(1:6),function(x){
    lines(s,r.n[,x],col=x,lwd=2)
})
legend("bottomleft",col=1:6,lwd=2,legend=c("core x core","far x core","salt x core","far x far","far x salt","salt x salt"))

# identify significant enrichments
# first remove FDR > 10%
for(i in 1:6){
    for(j in 1:length(s)){
        if(r.n[j,i+6] > 0.1){
            r.n[j,i] = NA
        }
    }
}
# then add a cross for significant ones
lapply(as.list(1:6),function(x){
    points(s,r.n[,x],col=x,pch=4)
})
dev.off()






#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################



# trans first
d = tr

# load pdb file

pdb = read.delim("000-additional_data/003-distances_trans.txt",header=F)
# convert positions to the ones used here
pdb$V1 = pdb$V1 - 161
pdb$V3 = pdb$V3 - 285

# keep only the positions of interest
pdb = pdb[pdb$V1 %in% 1:32 & pdb$V3 %in% 1:32,c(1,3,5)]
names(pdb) = c("pos1","pos2","dist")
pdb = pdb[order(pdb$pos1,pdb$pos2),]


# merge with epistasis scores
d2 = merge(d,pdb,by=c("pos1","pos2"))




# pretty much like for panel A
s = seq(0,0.5,0.01)
r.p = do.call("rbind",mclapply(as.list(s),function(e){
    
    r = do.call("rbind",lapply(as.list(c(3,5,10)),function(i){
        m = table(data.frame(factor(d2$dist < i,levels=c(T,F)),factor(d2$mepi > e & d2$hit,levels=c(T,F))))
        f = fisher.test(m)
        c(log(f$estimate),f$p.value)
    }))
    
    r[,2] = p.adjust(r[,2],method="fdr")
    as.vector(r)
}))
r.p[r.p == -Inf] = min(r.p[,1:3][is.finite(r.p[,1:3])])
r.p[r.p == Inf] = max(r.p[,1:3][is.finite(r.p[,1:3])])


r.n = do.call("rbind",mclapply(as.list(s),function(e){
    
    r = do.call("rbind",lapply(as.list(c(3,5,10)),function(i){
        m = table(data.frame(factor(d2$dist < i,levels=c(T,F)),factor(d2$mepi < -e & d2$hit,levels=c(T,F))))
        f = fisher.test(m)
        c(log(f$estimate),f$p.value)
    }))
    
    r[,2] = p.adjust(r[,2],method="fdr")
    as.vector(r)
}))
r.n[r.n == -Inf] = min(r.n[,1:3][is.finite(r.n[,1:3])])
r.n[r.n == Inf] = max(r.n[,1:3][is.finite(r.n[,1:3])])




pdf("Fig5S8C_trans.pdf")

plot(1,type="n",xlim=range(s),ylim=range(r.p[,1:3]),xlab="epistasis magnitude threshold",ylab="log10(odds ratio)")
lapply(as.list(1:3),function(x){
    lines(s,r.p[,x],col=x,lwd=2)
})
legend("bottomleft",col=1:3,lwd=2,legend=c("<3A","<5A","<10A"))
for(i in 1:3){
    for(j in 1:length(s)){
        if(r.p[j,i+3] > 0.1){
            r.p[j,i] = NA
        }
    }
}
lapply(as.list(1:3),function(x){
    points(s,r.p[,x],col=x,pch=4)
})


plot(1,type="n",xlim=range(s),ylim=range(r.n[,1:3]),xlab="epistasis magnitude threshold",ylab="log10(odds ratio)")
lapply(as.list(1:3),function(x){
    lines(s,r.n[,x],col=x,lwd=2)
})
legend("bottomleft",col=1:3,lwd=2,legend=c("<3A","<5A","<10A"))
for(i in 1:3){
    for(j in 1:length(s)){
        if(r.n[j,i+3] > 0.1){
            r.n[j,i] = NA
        }
    }
}
lapply(as.list(1:3),function(x){
    points(s,r.n[,x],col=x,pch=4)
})

dev.off()








# cis
d = ci

# load pdb file

pdb = read.delim("000-additional_data/003-distances_cis.txt",header=F)
# convert positions to the ones used here
pdb$V1 = pdb$V1 - 161
pdb$V3 = pdb$V3 - 161

# keep only the positions of interest
pdb = pdb[pdb$V1 %in% 1:32 & pdb$V3 %in% 1:32,c(1,3,5)]
names(pdb) = c("pos1","pos2","dist")
pdb = pdb[order(pdb$pos1,pdb$pos2),]


# merge with epistasis scores
d2 = merge(d,pdb,by=c("pos1","pos2"))




# pretty much like for panel A
s = seq(0,0.5,0.01)
r.p = do.call("rbind",mclapply(as.list(s),function(e){
    
    r = do.call("rbind",lapply(as.list(c(3,5,10)),function(i){
        m = table(data.frame(factor(d2$dist < i,levels=c(T,F)),factor(d2$mepi > e & d2$hit,levels=c(T,F))))
        f = fisher.test(m)
        c(log(f$estimate),f$p.value)
    }))
    
    r[,2] = p.adjust(r[,2],method="fdr")
    as.vector(r)
}))
r.p[r.p == -Inf] = min(r.p[,1:3][is.finite(r.p[,1:3])])
r.p[r.p == Inf] = max(r.p[,1:3][is.finite(r.p[,1:3])])


r.n = do.call("rbind",mclapply(as.list(s),function(e){
    
    r = do.call("rbind",lapply(as.list(c(3,5,10)),function(i){
        m = table(data.frame(factor(d2$dist < i,levels=c(T,F)),factor(d2$mepi < -e & d2$hit,levels=c(T,F))))
        f = fisher.test(m)
        c(log(f$estimate),f$p.value)
    }))
    
    r[,2] = p.adjust(r[,2],method="fdr")
    as.vector(r)
}))
r.n[r.n == -Inf] = min(r.n[,1:3][is.finite(r.n[,1:3])])
r.n[r.n == Inf] = max(r.n[,1:3][is.finite(r.n[,1:3])])




pdf("Fig5S8C_cis.pdf")

plot(1,type="n",xlim=range(s),ylim=range(r.p[,1:3]),xlab="epistasis magnitude threshold",ylab="log10(odds ratio)")
lapply(as.list(1:3),function(x){
    lines(s,r.p[,x],col=x,lwd=2)
})
legend("bottomleft",col=1:3,lwd=2,legend=c("<3A","<5A","<10A"))
for(i in 1:3){
    for(j in 1:length(s)){
        if(r.p[j,i+3] > 0.1){
            r.p[j,i] = NA
        }
    }
}
lapply(as.list(1:3),function(x){
    points(s,r.p[,x],col=x,pch=4)
})


plot(1,type="n",xlim=range(s),ylim=range(r.n[,1:3]),xlab="epistasis magnitude threshold",ylab="log10(odds ratio)")
lapply(as.list(1:3),function(x){
    lines(s,r.n[,x],col=x,lwd=2)
})
legend("bottomleft",col=1:3,lwd=2,legend=c("<3A","<5A","<10A"))
for(i in 1:3){
    for(j in 1:length(s)){
        if(r.n[j,i+3] > 0.1){
            r.n[j,i] = NA
        }
    }
}
lapply(as.list(1:3),function(x){
    points(s,r.n[,x],col=x,pch=4)
})

dev.off()




