library(parallel)
library(rgl)
source("../000-functions.R")

#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################

load("../Process_data/000-data/006-epistasis_thermo_model/trans.Rdata")

sx = seq(1.1,2,0.01)
sy = seq(1.1,2,0.01)
sb = seq(0.3,1.2,0.001)

dummy = matrix(Inf,ncol=length(sx),nrow=length(sy))

# load the whole array of parameter searched
# for each value of the B parameter
p = mclapply(1:901,function(j){
  load(paste0("../Process_data/000-data/004-param_search/",j,".Rdata"))
  dummy[!upper.tri(dummy)] = rss
  # the x and y parameters are intervertible in the model. parameters were y < x were not tested because redondant. So copy the lower half of the matrix to the upper half to make the parameter space symmetrical
  dummy[upper.tri(dummy)] = t(dummy)[upper.tri(dummy)]
  # compute variance explained
  matrix(1-dummy/sum((d$do_ppi1 - mean(d$do_ppi1))^2 + (d$do_ppi2 - mean(d$do_ppi2))^2 + (d$do_ppi3 - mean(d$do_ppi3))^2), ncol=length(sx))
},mc.cores=8)
p = simplify2array(p)

plot3d(1, type="n", xlim=range(sx), ylim=range(sy), zlim=range(p), xlab="", axes=F)
axes3d( edges=c("x-+", "y-+", "z") )
box3d()

# plot only one every 10 values of B, otherwise the plot is too heavy
for(i in seq(1,length(sb),10)){
  surface3d(sx ,sy , p[,,i], rainbow(length(sb)*1.5)[length(sb):1][i])
}

w = which(p == max(p), arr.ind=T)
points3d(sx[w[1,1]],sy[w[1,2]],max(p), size=10)

writeWebGL("Fig3S2A",width=1900,height=1000)

pdf("Fig3S2A_legend.pdf")
image(1:91,1,t(matrix(seq(1,length(sb),10),nrow=1)), col=rainbow(length(sb)*1.5)[length(sb):1][seq(1,length(sb),10)], axes=F, ylab="", main = paste0("At/ABwt=",sx[w[1,1]], ", Bt/ABwt=",sx[w[1,2]], ", b/ABwt=",sb[w[1,3]]))
axis(1, at = c(1,31,61,91), labels = c(0.3, 0.6, 0.9, 1.2))
dev.off()




#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################

# Load parameters from the Monte-Carlo cross validation
p = do.call("rbind",mclapply(1:1000,function(j){
    load(paste0("../Process_data/000-data/005-CV_thermo/",j,".Rdata"))
    out
},mc.cores=8))

pdf("Fig3S2B.pdf")
plot(p[,4], p[,5], xlab="proportion of variance explained train set", ylab="proportion of variance explained test set")
dev.off()



#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################


pdf("Fig3S2C.pdf")
pairs(d[,paste0("epi",1:3)], upper.panel = panel.cor.pearson, lower.panel = panel.smooth, diag.panel = panel.hist40)
dev.off()



#####################################################################################################################
### Panel D  ########################################################################################################
#####################################################################################################################

pdf("Fig3S2D.pdf")
hist(d$p, xlab="p-value", ylab="Counts")
dev.off()



#####################################################################################################################
### Panel E  ########################################################################################################
#####################################################################################################################

# plot only the 100 first points and then every 100 data points otherwise the plot is too heavy for illustrator
v = c(1:99,seq(100,nrow(d),100))

pdf("Fig3S2E.pdf")
plot(d$p[v], d$fdr[v]+d$fdrsd[v]*1.96, type="l", lwd = 2, col="grey80", xlab="p-value", ylab="FDR")
lines(d$p[v], d$fdr[v]+d$fdrsd[v]*1.96, col="grey80", lwd=2)
lines(d$p[v], d$fdr[v], lwd=2)
abline(0,1)
dev.off()



#####################################################################################################################
### Panel F  ########################################################################################################
#####################################################################################################################

p_thresh = min(d$p[!d$hit])

pdf("Fig3S2F.pdf")
plot(d$mepi, -log10(d$p), xlab="Genetic interaction score", ylab="-log10(p-value)", pch=".")
abline(h=-log10(p_thresh), lwd=2, col=2)
abline(v=c(-0.1,0.1), col=2, lwd=2)
text(-0.6,6,"negative genetic interactions", pos=4)
text(0.5,6,"positive genetic interactions", pos=2)
text(-0.6, -log10(p_thresh)-0.4, paste0("p < ", round(p_thresh,digits=6),"\nFDR < 20%"), col=2, pos=4)
dev.off()



#####################################################################################################################
### Panel G  ########################################################################################################
#####################################################################################################################


pos_hits = rev(cumsum(rev(table(cut(d$mepi[d$hit],seq(0,0.6,0.01))))))
neg_hits = rev(cumsum(table(cut(d$mepi[d$hit],seq(0,-0.6,-0.01), right=T))))

pdf("Fig3S2G.pdf")
plot(seq(0,0.59,0.01), neg_hits / pos_hits, type="l", lwd=2, xlab="magnitude threshold", ylab=" negative / positive interactions")
dev.off()