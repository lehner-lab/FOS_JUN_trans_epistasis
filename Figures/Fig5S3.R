source("../000-functions.R")
library(parallel)


#####################################################################################################################
### Panel A  ########################################################################################################
#####################################################################################################################

load("../Process_data/000-data/006-epistasis_thermo_model/trans_restricted.Rdata")
tr = d
load("../Process_data/000-data/006-epistasis_thermo_model/cis_restricted.Rdata")
ci = d

# split into matrix of singles and doubles
tr_s1 = as.matrix(tr[,paste0("s1_ppi",1:3)])
tr_s2 = as.matrix(tr[,paste0("s2_ppi",1:3)])
tr_do = as.matrix(tr[,paste0("do_ppi",1:3)])
ci_s1 = as.matrix(ci[,paste0("s1_ppi",1:3)])
ci_s2 = as.matrix(ci[,paste0("s2_ppi",1:3)])
ci_do = as.matrix(ci[,paste0("do_ppi",1:3)])


# fitted parameters
load("../Process_data/000-data/004-fitted_parameters.Rdata")


# compute variance explained by model fitted on both libs
# predictions with the thermo model fitted on trans
pr_tr_tr = thermo_model_pred(tr_s1,tr_s2, p$trans.r[1], p$trans.r[2], p$trans.r[3])
pr_tr_ci = thermo_model_pred(ci_s1,ci_s2, p$trans.r[1], p$trans.r[2], p$trans.r[3])
# predictions with the thermo model fitted on cis
pr_ci_tr = thermo_model_pred(tr_s1,tr_s2, p$cis.r[1], p$cis.r[2], p$cis.r[3])
pr_ci_ci = thermo_model_pred(ci_s1,ci_s2, p$cis.r[1], p$cis.r[2], p$cis.r[3])


# propotion of variance explained
var_tr_tr = 1 - colSums((tr_do-pr_tr_tr)^2) / colSums((tr_do-colMeans(tr_do))^2)
var_tr_ci = 1 - colSums((ci_do-pr_tr_ci)^2) / colSums((ci_do-colMeans(ci_do))^2)
var_ci_tr = 1 - colSums((tr_do-pr_ci_tr)^2) / colSums((tr_do-colMeans(tr_do))^2)
var_ci_ci = 1 - colSums((ci_do-pr_ci_ci)^2) / colSums((ci_do-colMeans(ci_do))^2)


pdf("Fig5S3A.pdf")
m = cbind(c(var_tr_tr,var_ci_tr),c(var_tr_ci,var_ci_ci))
barplot(m,beside=T, col = rep(c("grey40","grey80"),each=3), names.arg=c("trans","cis"),ylim=c(0,1), ylab="proportion of variance explained")
legend("topright", fill=c("grey40","grey80"), legend=c("fitted on trans","fitted on cis"))
dev.off()



#####################################################################################################################
### Panel B  ########################################################################################################
#####################################################################################################################

# load the results of the fitting on the matched libraries
p = do.call("rbind",mclapply(1:1000,function(j){
    load(paste0("../Process_data/000-data/008-fit_thermo_on_matched_libraries/",j,".Rdata"))
    out
},mc.cores=8))
colnames(p) = c("x_tr","y_tr","B_tr","var_exp_tr","var_exp_tr_on_ci","rmsd_tr","rmsd_tr_on_ci","x_ci","y_ci","B_ci","var_exp_ci","var_exp_ci_on_tr","rmsd_ci","rmsd_ci_on_tr")


xl = range(p[,c("var_exp_tr","var_exp_tr_on_ci","var_exp_ci","var_exp_ci_on_tr")])

pdf("Fig5S3B.pdf")
par(mfrow=c(2,2))
hist(p[,"var_exp_tr"], xlim=xl, main="", xlab="", ylab = "fitted on the trans library")
hist(p[,"var_exp_tr_on_ci"], xlim=xl, main="", xlab="", ylab="")
hist(p[,"var_exp_ci_on_tr"], xlim=xl, main="", xlab="variance explained in the trans library", ylab = "fitted on the cis library")
hist(p[,"var_exp_ci"], xlim=xl, main="", xlab="variance explained in the cis library", ylab="")
dev.off()





#####################################################################################################################
### Panel C  ########################################################################################################
#####################################################################################################################


# fitted parameters
load("../Process_data/000-data/004-fitted_parameters.Rdata")
x = p["x","cis.trans"]
y = p["y","cis.trans"]
B = p["B","cis.trans"]


prop = do.call("rbind",mclapply(1:1000,function(i){
    load(paste0("../Process_data/000-data/007-matched_libs/trans",i,".Rdata"))
    load(paste0("../Process_data/000-data/007-matched_libs/cis",i,".Rdata"))
    
    # split into matrix of singles and doubles
    tr_s1 = as.matrix(tr[,paste0("s1_ppi",1:3)])
    tr_s2 = as.matrix(tr[,paste0("s2_ppi",1:3)])
    tr_do = as.matrix(tr[,paste0("do_ppi",1:3)])
    ci_s1 = as.matrix(ci[,paste0("s1_ppi",1:3)])
    ci_s2 = as.matrix(ci[,paste0("s2_ppi",1:3)])
    ci_do = as.matrix(ci[,paste0("do_ppi",1:3)])
    
    # compute epistasis score for all replicates
    mult_epi_tr = tr_do - tr_s1 * tr_s2
    mult_epi_ci = ci_do - ci_s1 * ci_s2
    thermo_epi_tr = tr_do - thermo_model_pred(tr_s1,tr_s2,x,y,B)
    thermo_epi_ci = ci_do - thermo_model_pred(ci_s1,ci_s2,x,y,B)

    # variance explained the models
    var_mult_tr = 1-colSums(mult_epi_tr^2)/colSums((tr_do-colMeans(tr_do))^2)
    var_thermo_tr = 1-colSums(thermo_epi_tr^2)/colSums((tr_do-colMeans(tr_do))^2)
    var_mult_ci = 1-colSums(mult_epi_ci^2)/colSums((ci_do-colMeans(ci_do))^2)
    var_thermo_ci = 1-colSums(thermo_epi_ci^2)/colSums((ci_do-colMeans(ci_do))^2)
    
    # non-random variance
    nr_var_tr = cor(tr_do)^2
    nr_var_tr = mean(nr_var_tr[lower.tri(nr_var_tr)])
    nr_var_ci = cor(ci_do)^2
    nr_var_ci = mean(nr_var_ci[lower.tri(nr_var_ci)])
    
    # variance left to explain
    var_left_tr = nr_var_tr - var_mult_tr
    var_left_ci = nr_var_ci - var_mult_ci
    
    # proportion of non-random variance not explained by the multiplicztive but explained by thermodynamics
    prop_tr = (var_thermo_tr-var_mult_tr) / var_left_tr
    prop_ci = (var_thermo_ci-var_mult_ci) / var_left_ci
    
    c(prop_tr,prop_ci)
    
},mc.cores=8))

m = c(mean(prop[,1:3]),mean(prop[,4:6]))
s = c(sd(prop[,1:3]),sd(prop[,4:6])) / sqrt(3000) * 1.96

m = matrix(c(m,1-m),ncol=2)

# not necessary to plot confidence intervals because they are neglectable (0.0028 for trans and 0.00058 for cis)
pdf("Fig5S3C.pdf")
barplot(t(m), names.arg=c("trans","cis"), legend=c("explained by thermodynamics","residual genetic interactions"), ylab="proportion of non-random variance not explained by mutliplicativty")
dev.off()
