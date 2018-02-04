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


# which double mutants are hits? Logical vector used to construct the contingency table
hp = factor(tr$mepi > 0.1 & tr$hit, levels=c(T,F))
hn = factor(tr$mepi < -0.1 & tr$hit, levels=c(T,F))

r = do.call("rbind",lapply(as.list(1:32),function(x){
    do.call("rbind",mclapply(as.list(1:32),function(y){
        # which double mutants correspond to that position pair?
        pos = factor(tr$pos1 == x & tr$pos2 == y, levels=c(T,F))
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
pos = sort(unique(substr(tr$id1,1,2)))
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
r2.trans = rbind(tmp1,tmp2)

# this will be needed for Fig5S8
save(r2.trans, file="Fig5S6_data_forS8_trans.Rdata")



# now cis


# which double mutants are hits? Logical vector used to construct the contingency table
hp = factor(ci$mepi > 0.1 & ci$hit, levels=c(T,F))
hn = factor(ci$mepi < -0.1 & ci$hit, levels=c(T,F))

r = do.call("rbind",lapply(as.list(1:31),function(x){
    do.call("rbind",mclapply(as.list((x+1):32),function(y){
        # which double mutants correspond to that position pair?
        pos = factor(ci$pos1 == x & ci$pos2 == y, levels=c(T,F))
        # Fisher tests
        fp = fisher.test(table(pos,hp))
        fn = fisher.test(table(pos,hn))
        c(x,y,fp$estimate,fp$p.value,fn$estimate,fn$p.value)
    },mc.cores=8))
}))


# FDR
r = cbind(r,p.adjust(r[,4],method="fdr"))
r = cbind(r,p.adjust(r[,6],method="fdr"))

# add position identity
pos = sort(unique(substr(tr$id1,1,2)))
names(pos) = 1:32
r = as.data.frame(r)
r = cbind(pos[r[,1]],pos[r[,2]],r[,3:ncol(r)])

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
r2.cis = rbind(tmp1,tmp2)

# this will be needed for Fig5S8
save(r2.cis, file="Fig5S6_data_forS8_cis.Rdata")






pdf("Fig5S6A.pdf")
barplot(log10(r2.trans$V3),names.arg=paste(r2.trans$V1,r2.trans$V2),las=2,col=c(brewer.pal(9,"RdPu")[7],brewer.pal(9,"YlGn")[8])[r2.trans$V6],border=NA, ylab="log10(odds ratio")
barplot(log10(r2.cis$V3),names.arg=paste(r2.cis$V1,r2.cis$V2),las=2,col=c(brewer.pal(9,"RdPu")[7],brewer.pal(9,"YlGn")[8])[r2.cis$V6],border=NA, ylab="log10(odds ratio")
dev.off()



#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################

# for trans first
d = tr

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
max.prop = ceiling(max(c(pos, neg),na.rm=T)/10)*10
color = colorRampPalette(c("white",brewer.pal(9,"Set1")[2]))(10)
br = seq(0,max.prop,length.out=length(color)+1)


pdf("Fig5S6B_trans.pdf")
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

h1 = hist(neg[neg>0], breaks=100, plot=F)
h2 = hist(pos[pos>0], plot=F, breaks=h1$breaks)
xl = range(br)
yl = range(h1$counts, h2$counts)
hist(pos[pos>0], xlim=xl, ylim=yl, ylab="No. position pairs", xlab="Percentage of significant interactions per position pairs", breaks=h1$breaks)
hist(neg[neg>0], xlim=xl, ylim=yl, ylab="No. position pairs", xlab="Percentage of significant interactions per position pairs", breaks=h1$breaks)

dev.off()








# for cis
d = ci

# count percentage of hits by position pairs
perc = do.call("rbind",lapply(1:31,function(i){
    do.call("rbind",mclapply((i+1):32,function(j){
        x = d[d$pos1 == i & d$pos2 == j,]
        c(
            nrow(x[x$mepi > 0.1 & x$hit,]) / nrow(x) * 100,
            nrow(x[x$mepi < -0.1 & x$hit,]) / nrow(x) * 100
        )
        
    },mc.cores=8))
}))

# rows are Fos positions and cols are Jun positions
pos = matrix(ncol=32,nrow=32)
pos[lower.tri(pos)] = perc[,1]
pos[upper.tri(pos)] = t(pos)[upper.tri(pos)]
neg = matrix(ncol=32,nrow=32)
neg[lower.tri(neg)] = perc[,2]
neg[upper.tri(neg)] = t(neg)[upper.tri(neg)]


# same color scale as for trans color scale

pdf("Fig5S6B_cis.pdf")
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

