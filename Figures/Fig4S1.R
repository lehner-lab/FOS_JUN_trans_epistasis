library(RColorBrewer)
library(parallel)
library(aqfig)



# epistasis thresholds
s = seq(0,max(abs(d$mepi)),0.01)


# degreee of genetic interaction for each single mutants (i.e. number of interaction established)

# for Fos single mutants
d2 = split(d[d$hit,],d$id1[d$hit])
d.pos.f = do.call("rbind",mclapply(d2,function(x){
    rev(cumsum(rev(table(cut(x$mepi,s)))))
},mc.cores=8))

d.neg.f = do.call("rbind",mclapply(d2,function(x){
    rev(cumsum(rev(table(cut(-x$mepi,s)))))
},mc.cores=8))


# for Jun single mutants
d2 = split(d[d$hit,],d$id2[d$hit])
d.pos.j = do.call("rbind",mclapply(d2,function(x){
    rev(cumsum(rev(table(cut(x$mepi,s)))))
},mc.cores=8))

d.neg.j = do.call("rbind",mclapply(d2,function(x){
    rev(cumsum(rev(table(cut(-x$mepi,s)))))
},mc.cores=8))



# bin single mutants by degree
m = max(c(d.neg.j,d.neg.f,d.pos.j,d.pos.f))
br = seq(0,m+5,5)

d.neg.f = apply(d.neg.f,2,function(x){log10(table(cut(x,br,include.lowest=T)))})
d.neg.j = apply(d.neg.j,2,function(x){log10(table(cut(x,br,include.lowest=T)))})
d.pos.f = apply(d.pos.f,2,function(x){log10(table(cut(x,br,include.lowest=T)))})
d.pos.j = apply(d.pos.j,2,function(x){log10(table(cut(x,br,include.lowest=T)))})

# keep only until epistasis score threshold of 0.5
d.neg.f = d.neg.f[,1:50]
d.pos.f = d.pos.f[,1:50]
d.neg.j = d.neg.j[,1:50]
d.pos.j = d.pos.j[,1:50]


m.bin = ceiling(max(c(d.neg.j,d.neg.f,d.pos.j,d.pos.f),na.rm=T)*10)

pdf("Fig4S1.pdf")
par(mar=c(2,2,2,2),mfrow=c(2,2), oma=c(5,5,5,5))

image(1:nrow(d.pos.f),1:ncol(d.pos.f),d.pos.f,xlab="",ylab="",axes=F,col=heat.colors(m.bin),zlim=c(0,m.bin/10), main = "positive genetic interaction")
axis(2,at=seq(0,50,10),labels=seq(0,0.5,0.1),tick=F,cex.axis=0.7,las=2)
box()

image(1:nrow(d.neg.f),1:ncol(d.neg.f),d.neg.f,xlab="",ylab="",axes=F,col=heat.colors(m.bin),zlim=c(0,m.bin/10), main = "negative genetic interaction")
box()
mtext("Fos", side =4, line=2)

image(1:nrow(d.pos.j),1:ncol(d.pos.j),d.pos.j,xlab="",ylab="",axes=F,col=heat.colors(m.bin),zlim=c(0,m.bin/10))
axis(2,at=seq(0,50,10),labels=seq(0,0.5,0.1),tick=F,cex.axis=0.7,las=2)
axis(1,at=seq(par("usr")[1],par("usr")[2],length.out=length(br)),labels=br,tick=F,cex.axis=0.7,las=2)
box()

image(1:nrow(d.neg.j),1:ncol(d.neg.j),d.neg.j,xlab="",ylab="",axes=F,col=heat.colors(m.bin),zlim=c(0,m.bin/10))
axis(1,at=seq(par("usr")[1],par("usr")[2],length.out=length(br)),labels=br,tick=F,cex.axis=0.7,las=2)
box()
mtext("Jun", side =4, line=2)

mtext("No. of significant genetic interaction per single mutant", side = 1, line=2, outer=T)
mtext("genetic interaction magnitude threshold", side = 2, line=2, outer=T)

plot(1,type="n",xlab="",ylab="",axes=F)
vertical.image.legend(zlim=c(0,m.bin/10),col=heat.colors(m.bin))

dev.off()


