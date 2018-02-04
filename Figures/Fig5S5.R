library(parallel)
library(RColorBrewer)

load("../Process_data/000-data/006-epistasis_thermo_model/trans_restricted.Rdata")
tr = d
load("../Process_data/000-data/006-epistasis_thermo_model/cis_restricted.Rdata")
ci = d
load("../Process_data/000-data/006-epistasis_thermo_model/trans_restricted_filter100.Rdata")
trR = d
load("../Process_data/000-data/006-epistasis_thermo_model/cis_restricted_filter100.Rdata")
ciR = d


#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################

plot_enrich = function(d){
    
    # type of position
    d$type1[substr(d$id1,2,2) %in% c("a","d")] = "core"
    d$type1[substr(d$id1,2,2) %in% c("e","g")] = "salt b."
    d$type1[substr(d$id1,2,2) %in% c("b","c","f")] = "far"
    d$type2[substr(d$id2,2,2) %in% c("a","d")] = "core"
    d$type2[substr(d$id2,2,2) %in% c("e","g")] = "salt b."
    d$type2[substr(d$id2,2,2) %in% c("b","c","f")] = "far"
    
    # which double mutants are hits? Logical vector used to construct the contingency table
    hp = factor(d$mepi > 0.1 & d$hit, levels=c(T,F))
    hn = factor(d$mepi < -0.1 & d$hit, levels=c(T,F))
    
    # sets the order on the plot
    type_pairs = cbind(c("far","far","far","salt b.","salt b.","core"),c("far","salt b.","core","salt b.","core","core"))
    
    # compute enrichments
    r = t(apply(type_pairs,1,function(x){
        y = x[2]
        x = x[1]
        # which double mutants correspond to that position pair?
        pos = factor((d$type1 == x & d$type2 == y) | (d$type2 == x & d$type1 == y), levels=c(T,F))
        # Fisher tests
        fp = fisher.test(table(pos,hp))
        fn = fisher.test(table(pos,hn))
        c(x,y,fn$estimate,fn$p.value,fp$estimate,fp$p.value)
    }))
    
    r = as.data.frame(rbind(r[,1:4],r[,c(1:2,5:6)]))
    r[,3] = as.numeric(r[,3])
    r[,4] = as.numeric(r[,4])
    
    fdr = p.adjust(r[,4],method="fdr")
    sig = cut(fdr,c(0,0.001,0.01,0.1,1), labels=c("***","**","*","n.s."))
    
    b = barplot(log2(r[,3]), horiz=T, col=c(brewer.pal(9,"RdPu")[7],brewer.pal(9,"YlGn")[8])[rep(1:2,each=6)], xlab="log2(odds ratio)")
    labels = paste0(r$V1, " x ", r$V2, ", ", sig)
    lab.pos = sapply(log2(r[,3]),function(x){if(x > 0){2}else{4}})
    text(0, b, labels, pos=lab.pos)
}

pdf("Fig5S5A.pdf")
plot_enrich(tr)
plot_enrich(ci)
plot_enrich(trR)
plot_enrich(ciR)
dev.off()





#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################


# load distances in trans
dist = read.delim("000-additional_data/003-distances_trans.txt",header=F)
# convert positions to the ones used here
dist$V1 = dist$V1 - 161
dist$V3 = dist$V3 - 285

# keep only the positions of interest
dist = dist[dist$V1 %in% 1:32 & dist$V3 %in% 1:32,c(1,3,5)]
names(dist) = c("pos1","pos2","dist")
dist_trans = dist[order(dist$pos1,dist$pos2),]


# load distances in cis
dist = read.delim("000-additional_data/003-distances_cis.txt",header=F)
# convert positions to the ones used here
dist$V1 = dist$V1 - 161
dist$V3 = dist$V3 - 161

# keep only the positions of interest
dist = dist[dist$V1 %in% 1:32 & dist$V3 %in% 1:32,c(1,3,5)]
names(dist) = c("pos1","pos2","dist")
dist_cis = dist[order(dist$pos1,dist$pos2),]




enrich_dist = function(d, dist){
    # merge with epistasis scores
    d2 = merge(d,dist,by=c("pos1","pos2"))
    
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
    
    par(xpd=T)
    b = barplot(r, horiz=T, col=c(brewer.pal(9,"RdPu")[7],brewer.pal(9,"YlGn")[8])[rep(1:2,each=4)], xlab="log2(odds ratio)")
    labels = paste0(rep(c(">10A","<10A","<5A","<3A"),2), ", ", sig)
    lab.pos = sapply(r,function(x){if(x > 0){2}else{4}})
    text(0, b, labels, pos=lab.pos)
    if(flag.inf){
        text(min(r),b[r==min(r)],"-Inf",srt=90,pos=2)
    }
}

pdf("Fig5S5B.pdf")
enrich_dist(tr,dist_trans)
enrich_dist(trR,dist_trans)
enrich_dist(ci,dist_cis)
enrich_dist(ciR,dist_cis)
dev.off()
