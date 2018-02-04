library(parallel)


sx = seq(1.1,2,0.01)
sy = seq(1.1,2,0.01)
sb = seq(0.3,1.2,0.001)

dummy = matrix(Inf,ncol=length(sx),nrow=length(sy))

# for the trans library + the pooled cis-trans
p = do.call("cbind",lapply(1:4,function(i){
  # for each value of the B parameter
  p = mclapply(1:901,function(j){
    file = 901*(i-1) + j
    load(paste0("000-data/004-param_search/",file,".Rdata"))
    dummy[!upper.tri(dummy)] = rss
    dummy
  },mc.cores=8)
  p = simplify2array(p)
  w = which(p == min(p), arr.ind = T)
  c(sx[w[1,2]], sy[w[1,1]], sb[w[1,3]], min(p))
}))

p = as.data.frame(p)
names(p) = c("trans","trans.r","cis.r","cis.trans")
rownames(p) = c("x","y","B","rss")

save(p, file="000-data/004-fitted_parameters.Rdata")
