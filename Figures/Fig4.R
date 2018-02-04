library(parallel)
library(RColorBrewer)
library(aqfig)


#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################

load("../Process_data/000-data/006-epistasis_thermo_model/trans.Rdata")

# count percentage of hits by position pairs
d2 = split(d, list(d$pos1,d$pos2)) # sorted by Jun positions first
perc = do.call("rbind",mclapply(d2,function(x){
    c(
        nrow(x[x$mepi > 0.1 & x$hit,]) / nrow(x) * 100,
        nrow(x[x$mepi < -0.1 & x$hit,]) / nrow(x) * 100
    )
        
},mc.cores=8))

# rows are Fos positions and cols are Jun positions
pos = matrix(perc[,1],ncol=32)
neg = matrix(perc[,2],ncol=32)

# color scale
max.prop = ceiling(max(c(pos, neg))/10)*10
color = colorRampPalette(c("white",brewer.pal(9,"Set1")[2]))(10)
br = seq(0,max.prop,length.out=length(color)+1)


pdf("Fig4A.pdf")
par(mar=c(4,4,4,4))

image(1:32, 1:32, t(pos[32:1,]), col = color, breaks = br, axes=F, xlab="", ylab="")
axis(3,at=c(1,8,15,22,29,33)-0.5, tick=T, labels=NA)
axis(2,at=33-c(1,8,15,22,29,33)+0.5, tick=T, labels=NA)
box()
vertical.image.legend(zlim=range(br), col=color)

image(1:32, 1:32, t(neg[32:1,]), col = color, breaks = br, axes=F, xlab="", ylab="")
axis(3,at=c(1,8,15,22,29,33)-0.5, tick=T, labels=NA)
axis(2,at=33-c(1,8,15,22,29,33)+0.5, tick=T, labels=NA)
box()
vertical.image.legend(zlim=range(br), col=color)

h1 = hist(pos[pos>0], breaks=100, plot=F)
h2 = hist(neg[neg>0], plot=F, breaks=h1$breaks)
xl = range(br)
yl = range(h1$counts, h2$counts)
hist(pos[pos>0], xlim=xl, ylim=yl, ylab="No. position pairs", xlab="Percentage of significant interactions per position pairs", breaks=h1$breaks)
hist(neg[neg>0], xlim=xl, ylim=yl, ylab="No. position pairs", xlab="Percentage of significant interactions per position pairs", breaks=h1$breaks)

dev.off()




#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################


# which double mutants are hits? Logical vector used to construct the contingency table
hp = factor(d$mepi > 0.1 & d$hit, levels=c(T,F))
hn = factor(d$mepi < -0.1 & d$hit, levels=c(T,F))

r = do.call("rbind",lapply(as.list(1:32),function(x){
    do.call("rbind",mclapply(as.list(1:32),function(y){
        # which double mutants correspond to that position pair?
        pos = factor(d$pos1 == x & d$pos2 == y, levels=c(T,F))
        # Fisher tests
        fp = fisher.test(table(pos,hp))
        fn = fisher.test(table(pos,hn))
        c(fp$estimate,fp$p.value,fn$estimate,fn$p.value)
    },mc.cores=8))
}))


# FDR
r = cbind(r,p.adjust(r[,2],method="fdr"))
r = cbind(r,p.adjust(r[,4],method="fdr"))

# add position identity
pos = sort(unique(substr(d$id1,1,2)))
r = as.data.frame(r)
r = cbind(rep(pos,each=32),rep(pos,32),r)

# keep only significantly enriched position at FDR<10%
r.pos = r[r[,7] < 0.1,]
r.neg = r[r[,8] < 0.1,]
r.pos = r.pos[order(-r.pos[,3]),]
r.neg = r.neg[order(-r.neg[,5]),]

r.pos$color = 2
r.neg$color = 1

# combine the results for positive and negative
tmp1 = r.pos[,c(1:4,7,9)]
tmp2 = r.neg[,c(1:2,5:6,8:9)]
names(tmp1) = names(tmp2) = c("V1","V2","V3","V4","V5","V6")
r2 = rbind(tmp1,tmp2)

# this will be needed for Fig4S2
save(r2, file="Fig4_data_forS2.Rdata")

pdf("Fig4B.pdf")
barplot(log10(r2$V3),names.arg=paste(r2$V1,r2$V2),las=2,col=c(brewer.pal(9,"RdPu")[7],brewer.pal(9,"YlGn")[8])[r2$V6],border=NA, ylab="log10(odds ratio")
dev.off()



#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################

# type of position
d$type1[substr(d$id1,2,2) %in% c("a","d")] = "core"
d$type1[substr(d$id1,2,2) %in% c("e","g")] = "salt b."
d$type1[substr(d$id1,2,2) %in% c("b","c","f")] = "far"
d$type2[substr(d$id2,2,2) %in% c("a","d")] = "core"
d$type2[substr(d$id2,2,2) %in% c("e","g")] = "salt b."
d$type2[substr(d$id2,2,2) %in% c("b","c","f")] = "far"

types = unique(d$type1)


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
r = r[c(12,11,9,10,8,7,6,5,3,4,2,1),]

fdr = p.adjust(r[,4],method="fdr")
sig = cut(fdr,c(0,0.001,0.01,0.1,1), labels=c("***","**","*","n.s."))

pdf("Fig4C.pdf")
b = barplot(log2(r[,3]), horiz=T, col=c(brewer.pal(9,"RdPu")[7],brewer.pal(9,"YlGn")[8])[rep(1:2,each=6)], xlab="log2(odds ratio)")
labels = paste0(r$V1, " x ", r$V2, ", ", sig)
lab.pos = sapply(log2(r[,3]),function(x){if(x > 0){2}else{4}})
text(0, b, labels, pos=lab.pos)
dev.off()



#####################################################################################################################
### Panel D  ########################################################################################################
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
names(r) = NULL

# flag enrichemnts that are -Inf because there's no observations and set them to the minimal enrichment observed
if(length(r[!is.finite(r)])){
    flag.inf = T
}else{
    flag.inf = F
}
r[r == -Inf] = min(r[is.finite(r)])-0.1

# FDR and signifcance stars
fdr = p.adjust(c(fs.n$p.value,f10.n$p.value,f5.n$p.value,f3.n$p.value,fs$p.value,f10$p.value,f5$p.value,f3$p.value), method="fdr")
sig = cut(fdr,c(0,0.001,0.01,0.1,1), labels=c("***","**","*","n.s."))

pdf("Fig4D.pdf")
par(xpd=T)
b = barplot(r, horiz=T, col=c(brewer.pal(9,"RdPu")[7],brewer.pal(9,"YlGn")[8])[rep(1:2,each=4)], xlab="log2(odds ratio)")
labels = paste0(rep(c(">10A","<10A","<5A","<3A"),2), ", ", sig)
lab.pos = sapply(r,function(x){if(x > 0){2}else{4}})
text(0, b, labels, pos=lab.pos)
if(flag.inf){
    text(min(r),b[r==min(r)],"-Inf",srt=90,pos=2)
}
dev.off()



