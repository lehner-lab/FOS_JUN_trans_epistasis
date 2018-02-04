library(RColorBrewer)
library(aqfig)
library(ggplot2)
library(gridExtra)


aa = c("A","L","I","V","F","W","Y","H","S","T","Q","N","D","E","K","R","M","C","G","P")
wt_seq_F = strsplit("TDTLQAETDQLEDEKSALQTEIANLLKEKEKL",split="")[[1]]
wt_seq_J = strsplit("IARLEEKVKTLKAQNSELASTANMLREQVAQL",split="")[[1]]


#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################


load("../Process_data/000-data/002-PPIscores/trans_Q20_readCountFiltered.Rdata")

# single mutants
sf = d[d$pos2 == 0,]
sf = sf[,c("id1","pos1","wt1","mut1","mppi")]
sj = d[d$pos1 == 0,]
sj= sj[,c("id2","pos2","wt2","mut2","mppi")]
names(sj) = names(sf)

# make the matrices
m = merge(expand.grid(mut1=aa,pos1=1:32,stringsAsFactors=F),sf,all=T)
m = merge(m,data.frame(mut1=aa,ord=20:1))
m = m[order(m$pos1,m$ord),]
mf = matrix(m$mppi,ncol=32)

m = merge(expand.grid(mut1=aa,pos1=1:32,stringsAsFactors=F),sj,all=T)
m = merge(m,data.frame(mut1=aa,ord=20:1))
m = m[order(m$pos1,m$ord),]
mj = matrix(m$mppi,ncol=32)



# color scale
max.br = ceiling(max(c(mf,mj),na.rm=T)*10)/10
min.br = floor(min(c(mf,mj),na.rm=T)*10)/10
br = c(seq(min.br,1,0.025),seq(1.025,max.br,0.025))
colscale = c(colorRampPalette(c(brewer.pal(9,"RdBu")[9],"grey90"))(length(br[br<1])),colorRampPalette(c("grey90",brewer.pal(9,"RdBu")[1]))(length(br[br>1])))



pdf("Fig2A.pdf",width=7/20*32)
par(mar=c(5,5*32/20,5,5*32/20))
image(t(mf),col=colscale,breaks=br,axes=F)
vertical.image.legend(range(br),colscale)
text(seq(0,1,length.out=32),1+1/38,wt_seq_F,xpd=T,pos=3)
text(0-1/62,seq(0,1,length.out=20),aa[20:1],xpd=T,pos=2)

image(t(mj),col=colscale,breaks=br,axes=F)
vertical.image.legend(range(br),colscale)
text(seq(0,1,length.out=32),1+1/38,wt_seq_J,xpd=T,pos=3)
text(0-1/62,seq(0,1,length.out=20),aa[20:1],xpd=T,pos=2)
dev.off()





#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################


d = rbind(sf,sj)
d$prot = rep(c("Fos","Jun"),c(nrow(sf),nrow(sj)))
d$type[substr(d$id1,2,2) %in% c("a","d")] = "Hydrophobic Core"
d$type[substr(d$id1,2,2) %in% c("e","g")] = "Salt-Bridge Positions"
d$type[substr(d$id1,2,2) %in% c("b","c","f")] = "Far Side"


# test differences between position types
p1 = sprintf("%.1e", t.test(d$mppi[d$type == "Hydrophobic Core"], d$mppi[d$type == "Salt-Bridge Positions"])$p.value)
p2 = sprintf("%.1e", t.test(d$mppi[d$type == "Hydrophobic Core"], d$mppi[d$type == "Far Side"])$p.value)
p3 = sprintf("%.1e", t.test(d$mppi[d$type == "Far Side"], d$mppi[d$type == "Salt-Bridge Positions"])$p.value)

p = data.frame(core.salt = p1, core.far = p2, far.salt = p3)
row.names(p) = "p-value Welsh t-test"

pdf("Fig2B.pdf")
ggplot(d, aes(x=type, y=mppi, fill=prot)) +
  geom_violin(scale = "width") +
  theme_classic() +
  ylab("PPI scores") +
  scale_fill_manual(values=c("firebrick","deepskyblue"))

plot(1,type="n",xlab="",ylab="",axes=F)
grid.table(p)
dev.off()







#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################

# this code generates a pymol script that has to be run within pymol to color the structure

# average PPI score per position
mf = aggregate(sf$mppi, list(sf$pos1), mean)
mj = aggregate(sj$mppi, list(sj$pos1), mean)

# Color for protein residues not part of the domain of interest (as fraction of R G and B)
defaut_color <- "[0.25,0.25,0.25]" 

possible_AAs<- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
AA1to3 <- c(A="ALA", C="CYS", D="ASP", E="GLU", F="PHE", G="GLY", H="HIS", I="ILE", K="LYS", L="LEU", M="MET", N="ASN", P="PRO", Q="GLN", R="ARG", S="SER", T="THR", V="VAL", W="TRP", Y="TYR")

# position numbers in structure file
range_structure_F <- 162:(162+31) # Overlap between structure file and structure
pos_offset_F <- 161
range_structure_J <- 286:(286+31) # Overlap between structure file and
pos_offset_J <- 285


pymol_scores_tF <- mf$x
names(pymol_scores_tF) = 1:32 + pos_offset_F
pymol_scores_tF = pymol_scores_tF[names(pymol_scores_tF) %in% range_structure_F]
pymol_scores_tJ <- mj$x
names(pymol_scores_tJ) = 1:32 + pos_offset_J
pymol_scores_tJ = pymol_scores_tJ[names(pymol_scores_tJ) %in% range_structure_J]


# use the same colscale as in panel A
colStructure_tF <- colscale[cut(pymol_scores_tF, br)] # The max stuff is to have the scale centered on 0.
colStructure_tJ <- colscale[cut(pymol_scores_tJ, br)] # The max stuff is to have the scale centered on 0.



sink("Fig2C.pml")

# General commands
cat("fetch ", "1fos", "\n", sep="")
cat("show_as cartoon\n")
cat("util.cbc(selection='(all)',first_color=2,quiet=1,legacy=0,_self=cmd)\n")

# F and J chains as dots. Comment lines below to disable
cat("select chainE, chain E\n")
cat("show_as cartoon, chainE\n")
cat("show sticks, chainE\n")
cat("color yellow, chainE\n")

cat("select chainF, chain F\n")
cat("show_as cartoon, chainF\n")
cat("show sticks, chainF\n")
cat("color magenta, chainF\n")

cat("select chainG, chain G\n")
cat("hide everything, chainG\n")
cat("select chainH, chain H\n")
cat("hide everything, chainH\n")
cat("select chainC, chain C\n")
cat("hide everything, chainC\n")
cat("select chainD, chain D\n")
cat("hide everything, chainD\n")

# Now colouring each individual residue
for (i in 1:length(range_structure_F)) {      
  n <- paste("selF",i,sep="")    
  cat("select ", n, ", /", "1fos", "//E/", AA1to3[wt_seq_F[i]], "`", names(pymol_scores_tF)[i], "\n", sep="")
  cat("set_color newcolF" ,i, ", [", sep="") ;cat(round(col2rgb(colStructure_tF[i])[,1] / 255, 3), sep=","); cat("]\n", sep="")
  cat("color newcolF", i,", ", n ,"\n", sep="")
}
for (i in 1:length(range_structure_J)) {      
  n <- paste("selJ",i,sep="")    
  cat("select ", n, ", /", "1fos", "//F/", AA1to3[wt_seq_J[i]], "`", names(pymol_scores_tJ)[i], "\n", sep="")
  cat("set_color newcolJ" ,i, ", [", sep="") ;cat(round(col2rgb(colStructure_tJ[i])[,1] / 255, 3), sep=","); cat("]\n", sep="")
  cat("color newcolJ", i,", ", n ,"\n", sep="")
}

sink()




#####################################################################################################################
### Panel D  ########################################################################################################
#####################################################################################################################

load("../Process_data/000-data/002-PPIscores/trans_Q20_readCountFiltered.Rdata")

# merge Fos and Jun singles
sf = d[d$pos2 == 0,]
sf = sf[,c("id1","pos1","wt1","mut1",paste0("ppi",1:3),"mppi")]
sj = d[d$pos1 == 0,]
sj= sj[,c("id2","pos2","wt2","mut2",paste0("ppi",1:3),"mppi")]
names(sj) = names(sf)

d = merge(sf,sj,by="id1")
conf.int.x = apply(d[,paste0("ppi",1:3,".x")],1,sd) / sqrt(3) * 1.96
conf.int.y = apply(d[,paste0("ppi",1:3,".y")],1,sd) / sqrt(3) * 1.96

pdf("Fig2D.pdf")
plot(d$mppi.x,d$mppi.y,pch=16,cex=1.3,col="darkblue",xlab="Fos single mutants PPI score",ylab="Jun single mutants PPI score",xlim=range(c(d$mppi.x,d$mppi.y)),ylim=range(c(d$mppi.x,d$mppi.y)),main=paste("R =", round(cor(d$mppi.x,d$mppi.y),digits=3)))
abline(0,1,col="grey80",lty=2)
segments(d$mppi.x-conf.int.x,d$mppi.y,d$mppi.x+conf.int.x,d$mppi.y, col="darkblue")
segments(d$mppi.x,d$mppi.y-conf.int.y,d$mppi.x,d$mppi.y+conf.int.y, col="darkblue")
dev.off()


#####################################################################################################################
### Panel E  ########################################################################################################
#####################################################################################################################
# merge Fos and Jun single averages
mf = aggregate(sf$mppi, list(substr(sf$id1,1,2)),mean)
mj = aggregate(sj$mppi, list(substr(sj$id1,1,2)),mean)
d = merge(mf,mj,by="Group.1")

pdf("Fig2E.pdf")
plot(d$x.x,d$x.y,pch=16,type="n",col="darkblue",xlab="Fos single mutants PPI score,\nposition average",ylab="Jun single mutants PPI score,\nposition average",xlim=range(c(d$x.x,d$x.y)),ylim=range(c(d$x.x,d$x.y)),main=paste("R =", round(cor(d$x.x,d$x.y),digits=3)))
text(d$x.x,d$x.y,d$Group.1,cex=1.3)
abline(0,1,col="grey80",lty=2)
dev.off()
