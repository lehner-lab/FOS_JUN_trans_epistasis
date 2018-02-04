library(parallel)

load("../Process_data/000-data/006-epistasis_thermo_model/trans.Rdata")


#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################

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





pdf("Fig4S3A.pdf")

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
### Panel B  ########################################################################################################
#####################################################################################################################


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




pdf("Fig4S3B.pdf")

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




#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################

# load trans data filtered at 100 input reads
load("../Process_data/000-data/006-epistasis_thermo_model/trans_filter100.Rdata")

# type of position
d$type1[substr(d$id1,2,2) %in% c("a","d")] = "core"
d$type1[substr(d$id1,2,2) %in% c("e","g")] = "salt b."
d$type1[substr(d$id1,2,2) %in% c("b","c","f")] = "far"
d$type2[substr(d$id2,2,2) %in% c("a","d")] = "core"
d$type2[substr(d$id2,2,2) %in% c("e","g")] = "salt b."
d$type2[substr(d$id2,2,2) %in% c("b","c","f")] = "far"

types = unique(d$type1)

# which double mutants are hits? Logical vector used to construct the contingency table
hp = factor(d$mepi > 0.1 & d$hit, levels=c(T,F))
hn = factor(d$mepi < -0.1 & d$hit, levels=c(T,F))

r = do.call("rbind",lapply(1:length(types),function(x){
    do.call("rbind",lapply(x:length(types),function(y){
        # which double mutants correspond to that position pair?
        pos = factor((d$type1 == types[x] & d$type2 == types[y]) | (d$type2 == types[x] & d$type1 == types[y]), levels=c(T,F))
        # Fisher tests
        fp = fisher.test(table(pos,hp))
        fn = fisher.test(table(pos,hn))
        c(types[x],types[y],fp$estimate,fp$p.value,fn$estimate,fn$p.value)
    }))
}))

r = as.data.frame(rbind(r[,1:4],r[,c(1:2,5:6)]))
r[,3] = as.numeric(r[,3])
r[,4] = as.numeric(r[,4])

# reorder 
r = r[c(10,8,11,7,9,12,4,2,5,1,3,6),]

fdr = p.adjust(r[,4],method="fdr")
sig = cut(fdr,c(0,0.001,0.01,0.1,1), labels=c("***","**","*","n.s."))


pdf("Fig4S3C1.pdf")
b = barplot(log2(r[,3]), horiz=T, col=c(brewer.pal(9,"RdPu")[7],brewer.pal(9,"YlGn")[8])[rep(1:2,each=6)], xlab="log2(odds ratio)")
labels = paste0(r$V1, " x ", r$V2, ", ", sig)
lab.pos = sapply(log2(r[,3]),function(x){if(x > 0){2}else{4}})
text(0, b, labels, pos=lab.pos)
dev.off()





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



# compute enrichments at different distance thresholds
di = 3
t3 = table(data.frame(factor(d2$dist < di,levels=c(T,F)),factor(d2$mepi > 0.1 & d2$hit,levels=c(T,F))))
f3 = fisher.test(t3)

di = 5
t5 = table(data.frame(factor(d2$dist < di,levels=c(T,F)),factor(d2$mepi > 0.1 & d2$hit,levels=c(T,F))))
f5 = fisher.test(t5)

di = 10
t10 = table(data.frame(factor(d2$dist < di,levels=c(T,F)),factor(d2$mepi > 0.1 & d2$hit,levels=c(T,F))))
f10 = fisher.test(t10)

di = 10
ts = table(data.frame(factor(d2$dist >= di,levels=c(T,F)),factor(d2$mepi > 0.1 & d2$hit,levels=c(T,F))))
fs = fisher.test(ts)

di = 3
t3.n = table(data.frame(factor(d2$dist < di,levels=c(T,F)),factor(d2$mepi < -0.1 & d2$hit,levels=c(T,F))))
f3.n = fisher.test(t3.n)

di = 5
t5.n = table(data.frame(factor(d2$dist < di,levels=c(T,F)),factor(d2$mepi < -0.1 & d2$hit,levels=c(T,F))))
f5.n = fisher.test(t5.n)

di = 10
t10.n = table(data.frame(factor(d2$dist < di,levels=c(T,F)),factor(d2$mepi < -0.1 & d2$hit,levels=c(T,F))))
f10.n = fisher.test(t10.n)

di = 10
ts.n = table(data.frame(factor(d2$dist >= di,levels=c(T,F)),factor(d2$mepi < -0.1 & d2$hit,levels=c(T,F))))
fs.n = fisher.test(ts.n)

r = log2(c(fs.n$estimate,f10.n$estimate,f5.n$estimate,f3.n$estimate,fs$estimate,f10$estimate,f5$estimate,f3$estimate))
if(length(r[!is.finite(r)])){
    flag.inf = T
}else{
    flag.inf = F
}
r[r == -Inf] = min(r[is.finite(r)])-0.1
names(r) = NULL
fdr = p.adjust(c(fs.n$p.value,f10.n$p.value,f5.n$p.value,f3.n$p.value,fs$p.value,f10$p.value,f5$p.value,f3$p.value), method="fdr")
sig = cut(fdr,c(0,0.001,0.01,0.1,1), labels=c("***","**","*","n.s."))

pdf("Fig4S3C2.pdf")
par(xpd=T)
b = barplot(r, horiz=T, col=c(brewer.pal(9,"RdPu")[7],brewer.pal(9,"YlGn")[8])[rep(1:2,each=4)], xlab="log2(odds ratio)")
labels = paste0(rep(c(">10A","<10A","<5A","<3A"),2), ", ", sig)
lab.pos = sapply(r,function(x){if(x > 0){2}else{4}})
text(0, b, labels, pos=lab.pos)
if(flag.inf){
    text(min(r),b[r==min(r)],"-Inf",srt=90,pos=2)
}
dev.off()





