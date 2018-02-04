library(ggplot2)
library(gplots)
library(parallel)

#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################

load("../Process_data/000-data/006-epistasis_thermo_model/trans_restricted.Rdata")
tr = d
load("../Process_data/000-data/006-epistasis_thermo_model/cis_restricted.Rdata")
ci = d

# the cis library has to be made symmetrical, because singles can be matched to both Fos and Jun singles in the trans library
tmp = ci
tmp = tmp[,c(2,1,4,3,6,5,8,7,15:20,9:14,21:26,31:34,27:30,35:46)]
names(tmp) = names(ci)
ci2 = rbind(tmp, ci)

d = merge(tr,ci2,by=c("id1","id2"))
lim = range(d$do_mppi.x, d$do_mppi.y)
co = round(cor(d$do_mppi.x, d$do_mppi.y), digits=2)

pdf("Fig5B.pdf")
ggplot(d, aes(do_mppi.x, do_mppi.y)) +
    geom_hex(aes(fill=log2(..count..))) + 
    scale_x_continuous(limits = lim) + 
    scale_y_continuous(limits = lim) +
    xlab("trans double mutants PPI scores") + 
    ylab("cis double mutants PPI scores") + 
    annotate("text", x=0.45, y=1.1, label= paste0("R=",co)) +
    geom_abline(slope = 1, intercept = 0) + 
    theme_classic() +
    theme(legend.position = c(0.1, 0.8), panel.background = element_rect(color="black"))
dev.off()



#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################

lim = range(tr$do_mppi, ci$do_mppi)
lim = c(floor(lim[1]/0.04)*0.04, ceiling(lim[2]/0.04)*0.04)

# load matched libraries and bin them
s = seq(lim[1], lim[2], 0.04)
b = mclapply(1:1000, function(x){
    # trans lib
    load(paste0("../Process_data/000-data/007-matched_libs/trans", x, ".Rdata"))
    b1 = table(cut(tr$do_mppi,s))
    # cis lib
    load(paste0("../Process_data/000-data/007-matched_libs/cis", x, ".Rdata"))
    b2 = table(cut(ci$do_mppi,s))
    cbind(b1,b2)
},mc.cores=8)
b = simplify2array(b)

tr.all = b[,1,]
ci.all = b[,2,]

# frequencies
tr.all = tr.all / colSums(tr.all)
ci.all = ci.all / colSums(ci.all)

tr.m = rowMeans(tr.all)
ci.m = rowMeans(ci.all)
tr.conf = apply(tr.all,1,sd) / sqrt(ncol(tr)) * 1.96
ci.conf = apply(ci.all,1,sd) / sqrt(ncol(ci)) * 1.96


# bins for original libraries
tr2 = table(cut(tr$do_mppi,s))
ci2 = table(cut(ci$do_mppi,s))
tr2 = tr2 / sum(tr2)
ci2 = ci2 / sum(ci2)


yl = range(tr2, ci2, tr.m+tr.conf, ci.m+ci.conf)


pdf("Fig5C.pdf")
b = barplot2(rbind(tr.m, ci.m),beside=T,plot.ci=T,ci.l=rbind(tr.m-tr.conf,ci.m-ci.conf),ci.u=rbind(tr.m+tr.conf,ci.m+ci.conf),space=c(0,0),legend.text=c("trans","cis"), names.arg=NULL, ylim=yl, ylab="Density", xlab="PPI score")
axis(1,at=seq(0,50,10),labels=seq(0.2,1.2,0.2))
lines(colMeans(b), tr2, lwd=4, col=2)
lines(colMeans(b), ci2, lwd=4, col=7)
dev.off()





#####################################################################################################################
### Panel D  ########################################################################################################
#####################################################################################################################

# random variance is estimated as the average squared correlation coefficient between the three replicates PPI scores
nr_var_tr = cor(tr[,paste0("do_ppi",1:3)])^2
nr_var_tr = mean(nr_var_tr[lower.tri(nr_var_tr)])
nr_var_ci = cor(ci[,paste0("do_ppi",1:3)])^2
nr_var_ci = mean(nr_var_ci[lower.tri(nr_var_ci)])


# variance explained by the multiplicative model
rss = colSums((tr[,paste0("do_ppi",1:3)] - tr[,paste0("s1_ppi",1:3)] * tr[,paste0("s2_ppi",1:3)])^2)
tss = rowSums((t(tr[,paste0("do_ppi",1:3)]) - colMeans(tr[,paste0("do_ppi",1:3)]))^2)
mult_var_tr = 1 - rss / tss

rss = colSums((ci[,paste0("do_ppi",1:3)] - ci[,paste0("s1_ppi",1:3)] * ci[,paste0("s2_ppi",1:3)])^2)
tss = rowSums((t(ci[,paste0("do_ppi",1:3)]) - colMeans(ci[,paste0("do_ppi",1:3)]))^2)
mult_var_ci = 1 - rss / tss


# variance left to explain
var_left_tr = nr_var_tr - mult_var_tr
var_left_ci = nr_var_ci - mult_var_ci


# variance explained by thermodynamics
rss = colSums((tr[,paste0("epi",1:3)])^2)
tss = rowSums((t(tr[,paste0("do_ppi",1:3)]) - colMeans(tr[,paste0("do_ppi",1:3)]))^2)
thermo_var_tr = 1 - rss / tss
thermo_var_tr = thermo_var_tr - mult_var_tr

rss = colSums((ci[,paste0("epi",1:3)])^2)
tss = rowSums((t(ci[,paste0("do_ppi",1:3)]) - colMeans(ci[,paste0("do_ppi",1:3)]))^2)
thermo_var_ci = 1 - rss / tss
thermo_var_ci = thermo_var_ci - mult_var_ci


# propotion of non-random variance not explained by the multiplicative model that is explained by the thermodynamic model
prop_tr = thermo_var_tr / var_left_tr
prop_ci = thermo_var_ci / var_left_ci

av_prop_tr = mean(prop_tr)
sem_tr = sd(prop_tr) / sqrt(3)

av_prop_ci = mean(prop_ci)
sem_ci = sd(prop_ci) / sqrt(3)

pdf("Fig5D.pdf")
b = barplot2(matrix(c(av_prop_tr, 1-av_prop_tr,av_prop_ci, 1-av_prop_ci),ncol=2), col = c("grey40","grey80"), names.arg= c("trans","cis"), ylab="proportion of non-random variance left by the multiplicative model explained by thermodynamics")
segments(b, c(av_prop_tr-sem_tr,av_prop_ci-sem_ci), b, c(av_prop_tr+sem_tr,av_prop_ci+sem_ci))
dev.off()



#####################################################################################################################
### Panel E  ########################################################################################################
#####################################################################################################################

# p-value thresholds at FDR<20% in the trans library
p_tr = max(tr$p[tr$hit])


# corresponding FDR thresholds in both libraries
fdr_tr = min(tr$fdr[tr$p > p_tr])
fdr_ci = min(ci$fdr[ci$p > p_tr])


# proportion of total discoveries
tot_tr_pos = nrow(tr[tr$p <= p_tr & tr$mepi > 0,]) / nrow(tr)
tot_ci_pos = nrow(ci[ci$p <= p_tr & ci$mepi > 0,]) / nrow(ci)
tot_tr_neg = nrow(tr[tr$p <= p_tr & tr$mepi < 0,]) / nrow(tr)
tot_ci_neg = nrow(ci[ci$p <= p_tr & ci$mepi < 0,]) / nrow(ci)


# proportion of true discoveries
td_tr_pos = tot_tr_pos * (1-fdr_tr)
td_ci_pos = tot_ci_pos * (1-fdr_ci)
td_tr_neg = tot_tr_neg * (1-fdr_tr)
td_ci_neg = tot_ci_neg * (1-fdr_ci)


# matrix for plotting
m_pos = matrix(c(td_tr_pos, tot_tr_pos-td_tr_pos, td_ci_pos, tot_ci_pos-td_ci_pos),ncol=2) * 100
m_neg = matrix(c(td_tr_neg, tot_tr_neg-td_tr_neg, td_ci_neg, tot_ci_neg-td_ci_neg),ncol=2) * 100

yl = c(0,max(colSums(cbind(m_pos,m_neg))))
pdf("Fig5E.pdf")
b = barplot2(m_pos, col = c("grey40","grey80"), names.arg= c("trans","cis"), main="positive genetic interactions", ylab="Genetic interaction (%)", ylim=yl)
b = barplot2(m_neg, col = c("grey40","grey80"), names.arg= c("trans","cis"), main="negative genetic interactions", ylab="Genetic interaction (%)", ylim=yl)
dev.off()