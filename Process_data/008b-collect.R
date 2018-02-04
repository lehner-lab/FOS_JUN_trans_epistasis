library(parallel)

p = do.call("rbind",mclapply(1:1000,function(j){
  load(paste0("000-data/008-fit_thermo_on_matched_libraries/",j,".Rdata"))
  out
},mc.cores=8))
# get the median of each parameter and var. explained and RMSD
p = as.data.frame(apply(p,2,median))


names(p) = "cis.trans"
rownames(p) = c("x_tr","y_tr","B_tr","var_exp_tr","var_exp_tr_on_ci","rmsd_tr","rmsd_tr_on_ci","x_ci","y_ci","B_ci","var_exp_ci","var_exp_ci_on_tr","rmsd_ci","rmsd_ci_on_tr")

save(p, file="000-data/008-parameters_thermo_matched_libraries.Rdata")
