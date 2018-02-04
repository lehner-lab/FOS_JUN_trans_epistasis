source("../000-functions.R")
library(parallel)


# load PPI scores
load("../Process_data/000-data/002-PPIscores/trans_Q20_readCountFiltered.Rdata")

# single mutants
sf = d[d$pos2 == 0,]
sf = sf[,c("id1","pos1","wt1","mut1","mppi")]
sj = d[d$pos1 == 0,]
sj= sj[,c("id2","pos2","wt2","mut2","mppi")]
names(sj) = names(sf)



#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################


pdf("Fig2S1A.pdf")
x1 = hist(c(sf$mppi,sj$mppi), plot = F)
x2 = hist(sf$mppi,breaks=x1$breaks, plot = F)
x3 = hist(sj$mppi,breaks=x1$breaks, plot = F)
plot(1,type="n",xlim=range(x1$mids),ylim=range(c(x2$counts,x3$counts)),xlab="PPI scores (bin  middle)",ylab="Counts")
lines(x2$mids,x2$counts,lwd=4,col="firebrick")
lines(x3$mids,x3$counts,lwd=4,col="deepskyblue")
legend("topleft",lwd=4,col=c("firebrick","deepskyblue"),legend=c("Fos","Jun"))
dev.off()



#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################

# first need to compute if variant is significantly different from 0, for the whole library including double (because of FDR correction)
ppi = as.matrix(d[,paste0("ppi",1:3)])
m = rowMeans(ppi)
s = sqrt(rowSums((ppi-m)^2)/2) / sqrt(3)
d$p = one.sample.ttest(m,s,1,3)
d$fdr = p.adjust(d$p,method="fdr")

# get singles
sf = d[d$pos2 == 0,]
sf = sf[,c("id1","pos1","wt1","mut1","mppi","p","fdr")]
sj = d[d$pos1 == 0,]
sj= sj[,c("id2","pos2","wt2","mut2","mppi","p","fdr")]
names(sj) = names(sf)

sf2 = sf[sf$fdr < 0.05,]
sj2 = sj[sj$fdr < 0.05,]

# proportion of detrimental variants at various thresholds
s = c(0,seq(0.4,0.8,0.01))
pf = cumsum(table(cut(sf2$mppi,s))) / nrow(sf)
pj = cumsum(table(cut(sj2$mppi,s))) / nrow(sj)


pdf("Fig2S1B.pdf")
plot(s[2:length(s)],pf,type="l",xlab="PPIscore threshold",ylab="proportion detrimental",col="firebrick",lwd=3,ylim=range(c(pf,pj)))
lines(s[2:length(s)],pj,col="deepskyblue",lwd=3)
legend("topleft",col=c("firebrick","deepskyblue"),lwd=3,legend=c("single Fos","single Jun"))
dev.off()



#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################

s = seq(1,1.1,0.001)
pf = rev(cumsum(rev(table(cut(sf2$mppi,s))))) / nrow(sf)
pj = rev(cumsum(rev(table(cut(sj2$mppi,s))))) / nrow(sj)

pdf("Fig2S1C.pdf")
plot(s[2:length(s)],pf,type="l",xlab="PPIscore threshold",ylab="proportion stronger",col="firebrick",lwd=3,ylim=range(c(pf,pj)))
lines(s[2:length(s)],pj,col="deepskyblue",lwd=3)
legend("topright",col=c("firebrick","deepskyblue"),lwd=3,legend=c("single Fos","single Jun"))
dev.off()




#####################################################################################################################
### Panel D  ########################################################################################################
#####################################################################################################################

# load amino acid features data
aa_indexes = read.delim("000-additional_data/002-index1_formatted.txt", header=F)
aa = strsplit("ARNDCQEGHILKMFPSTWYV",split="")[[1]]
names(aa_indexes) = c("id","index",aa)

# remove features with NAs
aa_indexes = aa_indexes[!is.na(rowSums(aa_indexes[,3:22])),]

# compute the change in the magnitude of each feature for each substitution (mut minus wt)
delta = do.call("cbind",mclapply(1:nrow(aa_indexes),function(x){
  x = as.numeric(as.character(aa_indexes[x,3:22]))
  # matrix of changes in magnitudes. Rows are mutants aa and columns are wt aa
  delta = outer(x,x,"-")
  
  as.vector(delta)
},mc.cores=8))
delta = data.frame(expand.grid(aa,aa, stringsAsFactors = F), delta)
names(delta) = c("mut1","wt1",aa_indexes[,1])


# merge with PPI scores
sf = merge(sf, delta, by=c("wt1","mut1"))
sj = merge(sj, delta, by=c("wt1","mut1"))


# linear regression by position
sf2 = split(sf, sf$pos1)
sj2 = split(sj, sj$pos1)

var_sf = do.call("rbind",lapply(sf2,function(y){
  do.call("c",mclapply(as.list(8:ncol(y)),function(x){
    summary(lm(y$mppi ~ y[,x]))$r.squared
  },mc.cores=8))
}))

var_sj = do.call("rbind",lapply(sj2,function(y){
  do.call("c",mclapply(as.list(8:ncol(y)),function(x){
    summary(lm(y$mppi ~ y[,x]))$r.squared
  },mc.cores=8))
}))


# max variance explained by any feature at each position
max_f = apply(var_sf,1,max)
max_j = apply(var_sj,1,max)


# order by heptad position (a-g)
ord = order(rep(1:7,5)[1:32]) # vector for reordering a vector intially sorted by position (1 - 32)
max_f = max_f[ord]
max_j = max_j[ord]
names(max_f) = sort(unique(paste0(substr(sf$id1,1,2),sf$wt1)))[ord]
names(max_j) = sort(unique(paste0(substr(sj$id1,1,2),sj$wt1)))[ord]


pdf("Fig2S2D.pdf")
par(mfrow=c(2,1))
barplot(max_f,las=2,col=c(rep(1:4,each=5),rep(5:7,each=4)),ylab="max var explained",main="Fos",ylim=c(0,1))
barplot(max_j,las=2,col=c(rep(1:4,each=5),rep(5:7,each=4)),ylab="max var explained",main="Jun",ylim=c(0,1))
dev.off()

