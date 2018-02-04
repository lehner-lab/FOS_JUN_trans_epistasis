library(parallel)

# for the trans library + the pooled cis-trans
p = do.call("cbind",lapply(1:4,function(i){
  # for each value of the B parameter
  p = do.call("rbind",mclapply(1:1000,function(j){
    file = 1000*(i-1) + j
    load(paste0("000-data/005-CV_thermo/",file,".Rdata"))
    out
  },mc.cores=8))
  # get the median of each parameter and var. explained and RMSD
  apply(p,2,median)
}))

p = as.data.frame(p)
names(p) = c("trans","trans.r","cis.r","cis.trans")
rownames(p) = c("x","y","B","var_exp_train","var_exp_test","rmsd")

save(p, file="000-data/005-crossValidated_parameters.Rdata")
