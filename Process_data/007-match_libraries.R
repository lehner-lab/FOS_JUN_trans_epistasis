
job = as.numeric(as.character(commandArgs(trailingOnly = TRUE)))

set.seed(job)

# number of bins into which the libraries are split
nbin = 20


# load data
load("000-data/003-mult_epi/trans_restricted.Rdata")
tr = d
load("000-data/003-mult_epi/cis_restricted.Rdata")
ci = d


# define bins
r = range(c(tr$s1_mppi,tr$s2_mppi,ci$s1_mppi,ci$s2_mppi))
s = seq(r[1],r[2],length.out=nbin+1)



# the trans library can already be binned
# for the cis library, single mutants need to be divided randomly within one or the other dimension
tr$bin1 = cut(tr$s1_mppi,s,include.lowest=T)
tr$bin2 = cut(tr$s2_mppi,s,include.lowest=T)
tr2 = split(tr,as.list(tr[,c("bin1","bin2")]))

# counts in each bin
counts_trans = do.call("c", lapply(tr2,nrow))

# 10 randomization for that job
sapply(1:10,function(i){
  
  # randomly assign cis single mutants to one or the other dimension
  rand = sample(c(T,F),nrow(ci),replace=T)
  tmp1 = ci[rand,]
  tmp2 = ci[!rand,]
  tmp2 = tmp2[,c(2,1,4,3,6,5,8,7,15:20,9:14,21:26,31:34,27:30,35:46)]
  names(tmp2) = names(tmp1)
  ci = rbind(tmp1,tmp2)
  
  # bin the cis library
  ci$bin1 = cut(ci$s1_mppi,s,include.lowest=T)
  ci$bin2 = cut(ci$s2_mppi,s,include.lowest=T)
  ci2 = split(ci,as.list(ci[,c("bin1","bin2")]))
  
  # counts in each bin
  counts_cis = do.call("c", lapply(ci2,nrow))
  
  # count differences between the two libraries
  counts_diff = counts_trans - counts_cis
  
  # for each bin, if there is a higher number of trans, sub-sample, and vice-versa
  for(j in seq_along(counts_diff)){
    if(counts_diff[j] > 0){
        tr2[[j]] = tr2[[j]][sample(counts_trans[j],counts_cis[j]),]
    }else if(counts_diff[j] < 0){
        ci2[[j]] = ci2[[j]][sample(counts_cis[j],counts_trans[j]),]
    }
  }
  
  tr = do.call("rbind",tr2)[,1:38]
  ci = do.call("rbind",ci2)[,1:38]
  
  out = (job-1) * 10 + i
  
  save(tr,file=paste0("000-data/007-matched_libs/trans",out,".Rdata"))
  save(ci,file=paste0("000-data/007-matched_libs/cis",out,".Rdata"))
  
})
