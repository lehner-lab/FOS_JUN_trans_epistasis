library(ggplot2)
library(RColorBrewer)
library(rgl)
source("../000-functions.R")

#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################

# load data
load("../Process_data/000-data/002-PPIscores/trans_Q20_readCountFiltered.Rdata")

# single mutants
sf = d[d$pos2 == 0,]
sf = sf[,c("id1","pos1","wt1","mut1","ppi1")]
sj = d[d$pos1 == 0,]
sj= sj[,c("id2","pos2","wt2","mut2","ppi1")]
names(sj) = names(sf)


# doubles
do = d[d$pos1 != 0 & d$pos2 != 0,]
# since there is one single mutant of Fos missing, remove doubles that contain it
do = do[do$id1 %in% sf$id1,]


pdf("Fig3A.pdf")
x=hist(do$ppi1,breaks=40,plot=F)
y=hist(sf$ppi1,plot=F,breaks=x$breaks)
y2=hist(sj$ppi1,plot=F,breaks=x$breaks)
plot(x,ylim=range(c(x$density,y$density,y2$density)),freq=F,border=NA,col="grey80",xlab="PPI score")
lines(y$mids,y$density,col="firebrick",lwd=5)
lines(y2$mids,y2$density,col="deepskyblue",lwd=5)
abline(v=1,lty=2)
legend("topleft",fill=c("grey80","firebrick","deepskyblue"),legend=c("double mutants","FOS single mutants","JUN single mutants"))
dev.off()



#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################

# load multiplicative epistasis data
load("../Process_data/000-data/003-mult_epi/trans.Rdata")

# temporary data.frame for plotting
tmp = data.frame(expected = d$s1_ppi1 * d$s2_ppi1, observed = d$do_ppi1)


pdf("Fig3B.pdf")
ggplot(tmp, aes(expected,observed)) +
  geom_hex(aes(fill=log10(..count..))) + 
  scale_x_continuous(limits = range(tmp)) + 
  scale_y_continuous(limits = range(tmp)) +
  geom_abline(slope = 1, intercept = 0) + 
  theme_classic() +
  theme(legend.position = c(0.1, 0.8), panel.background = element_rect(color="black"))
dev.off()



#####################################################################################################################
### Panel C and F ###################################################################################################
#####################################################################################################################

# function for the plot
pie_chart = function(d, lib){

  # bin epistasis scores
  d$bin.epi = cut(d$epi1, br.epi, include.lowest=T)
  
  # split the data into bins of single mutants
  d$bin1 = cut(d$s1_ppi1, br.singles, include.lowest=T)
  d$bin2 = cut(d$s2_ppi1, br.singles, include.lowest=T)
  d2 = split(d,list(d$bin1,d$bin2))

  # layout of the plot. Will be plotting by row
  mat = matrix(1:length(d2),ncol=length(br.singles)-1,byrow=T)
  # the rows have to be inverted because we start plotting the bins with the smallest PPI scores
  mat = mat[nrow(mat):1,]
  # add some space around for labels and legend
  mat = cbind(length(d2)+1,mat,length(d2)+3,length(d2)+3)
  mat = rbind(length(d2)+4,length(d2)+4,mat,c(0,rep(length(d2)+2,ncol(mat)-3),0,0))
  par(oma=c(0,0,0,0),mar=c(0,0,0,0),xpd=T)
  layout(mat)
  
  # for each two-way bin
  for(i in 1:length(d2)){
    if(nrow(d2[[i]]) > 0){
      # compute the proportion of each epistasis bin and then plot the pie chart
      prop = table(d2[[i]]$bin.epi)
      pie(prop,col=colscale,labels=NA,clockwise=T,border=NA,radius=1)         
    }else{
      plot(1,type="n",xlab="",ylab="",axes=F)
    }
  }
  
  # axis abels
  par(mar=c(0.5,0,0.5,0))
  plot(1,type="n",xlab="",ylab="",axes=F)
  text(1.3,seq(par("usr")[3],par("usr")[4],length.out=length(br.singles)),labels=br.singles)
  text(0.8,1,labels="PPI score Jun single mutant",srt=90,cex=2)
  
  par(mar=c(0,0.5,0,0.5))
  plot(1,type="n",xlab="",ylab="",axes=F)
  text(seq(par("usr")[3],par("usr")[4],length.out=length(br.singles)),1.3,labels=br.singles)
  text(1,0.8,labels="PPI score Fos single mutant",cex=2)   
  
  # legend
  par(mar=c(1,2,1,3))
  image(matrix(1:length(colscale),nrow=1),col=colscale,axes=F)
  axis(4,at=seq(par("usr")[3],par("usr")[4],length.out=length(br.epi)),labels=round(br.epi,digits=1),las=2)
  mtext("epistasis score",side=2)
  
  # main title
  plot(1,type="n",xlab="",ylab="",axes=F)
  text(1,1,labels=lib,cex=2)
  
}



# epistasis data
mult = d
load("../Process_data/000-data/006-epistasis_thermo_model/trans.Rdata")
thermo = d

# single mutants bins
max.ppi = ceiling(max(c(mult$s1_ppi1,mult$s2_ppi1,thermo$s1_ppi1,thermo$s2_ppi1))*10) / 10
min.ppi = floor(min(c(mult$s1_ppi1,mult$s2_ppi1,thermo$s1_ppi1,thermo$s2_ppi1))*10) / 10
br.singles = seq(min.ppi,max.ppi,0.05)

# epistasis bins
max.epi = ceiling(max(abs(c(mult$epi1,thermo$epi1)))*10) / 10
br.epi = seq(-max.epi,max.epi,length.out=max.epi*20+1)

# color scale for epistasis bins
colscale = colorRampPalette(c(brewer.pal(9,"RdPu")[7],"#FFFFBF",brewer.pal(9,"YlGn")[8]))(max.epi*20)


pdf("Fig3C_F.pdf")
pie_chart(mult, "multiplicative model")
pie_chart(thermo, "thermodynamic model")
dev.off()



#####################################################################################################################
### Panel D  ########################################################################################################
#####################################################################################################################

# load thermodynamic model parameters
load("../Process_data/000-data/004-fitted_parameters.Rdata")
x = p$trans[1]
y = p$trans[2]
B = p$trans[3]

min.ppi = B/(B+1)
max.ppi = (min(c(x,y)) + B)/ (B + 1)

# function to comput ddG from ppi
pred_ppi = function(ddG,x,y,B){
    X = (B+x)/(B+1)
    Y = (B+y)/(B+1)
    
    b = ((X-1)*(Y-1)*(B+1)*exp(ddG)) + X + Y
    
    0.5 * (b - sqrt(b^2-4*(B*(X-1)*(Y-1)*exp(ddG)+X*Y)))
}

s = seq(-10,15,0.1)

pdf("Fig3D.pdf")
plot(s,pred_ppi(s,x,y,B), type="l", xlab="ddG (a.u.)", ylab="PPI score")
dev.off()





#####################################################################################################################
### Panel E  ########################################################################################################
#####################################################################################################################


plot3d(   x = thermo$s1_ppi1,
          y = thermo$s2_ppi1,
          z = thermo$do_ppi1,
          xlab = "", ylab = "", zlab = "",
          col = colscale[cut(thermo$epi1, br.epi)],
          xlim= range(c(thermo$s1_ppi1,thermo$s2_ppi1)),
          ylim= range(c(thermo$s1_ppi1,thermo$s2_ppi1))
)


ranges <- rgl:::.getRanges()
s1 <- seq(ranges$xlim[1], ranges$xlim[2], length=30)
s2 <- seq(ranges$xlim[1], ranges$xlim[2], length=30)
do <- outer(s1,s2,function(s1,s2){thermo_model_pred(s1,s2,x,y,B)})


surface3d(s1, s2, do, alpha=0.5, lit=F)

writeWebGL("Fig3E",width=1900,height=1000)




#####################################################################################################################
### Panel F  ########################################################################################################
#####################################################################################################################

# ploted with 3C



#####################################################################################################################
### Panel G  ########################################################################################################
#####################################################################################################################

# compute average epistasis score by single mutants
mf = aggregate(mult$epi1, as.list(mult[,c("id1","s1_ppi1")]), mean)
mj = aggregate(mult$epi1, as.list(mult[,c("id2","s2_ppi1")]), mean)
tf = aggregate(thermo$epi1, as.list(thermo[,c("id1","s1_ppi1")]), mean)
tj = aggregate(thermo$epi1, as.list(thermo[,c("id2","s2_ppi1")]), mean)

mf = mf[order(mf$s1_ppi1),]
mj = mj[order(mj$s2_ppi1),]
tf = tf[order(tf$s1_ppi1),]
tj = tj[order(tj$s2_ppi1),]


xl = range(mf$s1_ppi1, mj$s2_ppi1)
yl = range(mf$x, mj$x, tf$x, tj$x)



pdf("Fig3G.pdf",width=14)
par(mfrow=c(2,2), mar=c(0,0,0,0), oma=c(5,5,5,5), xpd=T)

plot(mf$s1_ppi1, mf$x, pch=mf$id1, axes=F, xlab="", ylab="", xlim = xl, ylim = yl, type = "n")
text(mf$s1_ppi1, mf$x, mf$id1)
abline(h=0, col=2, lwd=2)
lines(mf$s1_ppi1, loess(mf$x ~ mf$s1_ppi1)$fitted, col=5, lwd=2)
axis(2)
mtext("Fos", side=3)
box()

plot(mj$s2_ppi1, mj$x, pch=mj$id2, axes=F, xlab="", ylab="", xlim = xl, ylim = yl, type = "n")
text(mj$s2_ppi1, mj$x, mj$id2)
abline(h=0, col=2, lwd=2)
lines(mj$s2_ppi1, loess(mj$x ~ mj$s2_ppi1)$fitted, col=5, lwd=2)
mtext("multiplicative", side=4, line=2)
mtext("Jun", side=3)
box()

plot(tf$s1_ppi1, tf$x, pch=tf$id1, xlab="", ylab="", xlim = xl, ylim = yl, type = "n")
text(tf$s1_ppi1, tf$x, tf$id1)
abline(h=0, col=2, lwd=2)
lines(tf$s1_ppi1, loess(tf$x ~ tf$s1_ppi1)$fitted, col=5, lwd=2)

plot(tj$s2_ppi1, tj$x, pch=tj$id2, axes=F, xlab="", ylab="", xlim = xl, ylim = yl, type = "n")
text(tj$s2_ppi1, tj$x, tj$id2)
abline(h=0, col=2, lwd=2)
lines(tj$s2_ppi1, loess(tj$x ~ tj$s2_ppi1)$fitted, col=5, lwd=2)
axis(1)
box()

mtext("thermodynamic", side=4, line=2)
mtext("single mutant PPI score", side=1, line=3, outer=T)
mtext("Average genetic interaction score", side=2, line=3, outer=T)
dev.off()




