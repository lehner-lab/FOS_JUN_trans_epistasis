job = as.numeric(as.character(commandArgs(trailingOnly = TRUE)))
set.seed(job)

#####################################################################################################################
### Compute epistasis scores ########################################################################################
#####################################################################################################################


source("../000-functions.R")

# load data
load(paste0("000-data/007-matched_libs/trans",job,".Rdata"))
load(paste0("000-data/007-matched_libs/cis",job,".Rdata"))

# get thermodynamic model parameters
load("000-data/004-fitted_parameters.Rdata")
x = p["x","cis.trans"]
y = p["y","cis.trans"]
B = p["B","cis.trans"]



# compute epistasis score for all replicates
epi_trans = sapply(1:3,function(i){
  pr = thermo_model_pred(tr[,paste0("s1_ppi",i)],tr[,paste0("s2_ppi",i)],x,y,B)
  tr[,paste0("do_ppi",i)] - pr
})
epi_trans = cbind(epi_trans,rowMeans(epi_trans))
stderr = sqrt(rowSums((epi_trans[,1:3]-epi_trans[,4])^2)/2) / sqrt(3)
epi_trans = cbind(epi_trans, one.sample.ttest(epi_trans[,4],stderr,0,3))

# same for the cis library
epi_cis = sapply(1:3,function(i){
  pr = thermo_model_pred(ci[,paste0("s1_ppi",i)],ci[,paste0("s2_ppi",i)],x,y,B)
  ci[,paste0("do_ppi",i)] - pr
})
epi_cis = cbind(epi_cis,rowMeans(epi_cis))
stderr = sqrt(rowSums((epi_cis[,1:3]-epi_cis[,4])^2)/2) / sqrt(3)
epi_cis = cbind(epi_cis, one.sample.ttest(epi_cis[,4],stderr,0,3))



#####################################################################################################################
### Compute FDR by permutations #####################################################################################
#####################################################################################################################


permut = function(d, out_file){
  
  # order by p-values
  d = d[order(d[,5]),]
  nd = nrow(d)
  rp = d[,5]
  
  # for each randomization
  # (need to split it in chunks of 100 because of memory issues)
  q = do.call("cbind",lapply(1:(nPermut/100),function(x){
    cat("\r", x, "\t\t\t")
    q = sapply(1:100,function(y){
      
      # randomize
      d[,1:3] = apply(d[,1:3],2,sample)
      
      # compute mean and standard error
      d[,4] = rowMeans(d[,1:3])
      s = sqrt(rowSums((d[,1:3]-d[,4])^2)/2) / sqrt(3)
      
      # t-test
      p = one.sample.ttest(d[,4],s,0,3)
      
      # count number of false discoveries at each p-value threshold
      q = cut(p,breaks=c(0,rp))
      
      # cumulative sum to count all the false discoveries below a given p-value threshold
      cumsum(table(q))
      
    })
    
    # average number of false discoveries and standard error at each p-value threshold across this batch of 100 permutations
    m = rowMeans(q)
    s = sqrt(rowSums((q-m)^2)/(nPermut-1) / nPermut)
    
    cbind(m,s)
    
  }))
  
  # average number of false discoveries and standard error at each p-value threshold across all permutations
  m = rowMeans(q[,seq(1,nPermut/100*2,2)])
  s = sqrt(rowSums(q[,seq(2,nPermut/100*2,2)]^2))
  
  # number of total discoveries at each p-value threshold
  o = 1:length(rp)
  
  # false discovery rate
  fdr = m/o
  
  
  # total discoveries at this p-value threshold
  # for positive and negative epistasis separately
  td_pos = cumsum(table(cut(d[d[,4] > 0,5],c(0,seq_p))))
  td_neg = cumsum(table(cut(d[d[,4] < 0,5],c(0,seq_p))))
  
  # fdr at all p-value threshold
  f = sapply(seq_p,function(x){
    max(fdr[d[,5] < x])
  })
  
  # proportion of true discoveries
  prop_pos = td_pos * (1-f) / nrow(d)
  prop_neg = td_neg * (1-f) / nrow(d)
  prop = cbind(prop_pos,prop_neg)

  save(prop,file=out_file)
}



# sequence of p-value thresholds at which to compute the proportion of true discoveries
seq_p = 10^seq(-5,0,length.out=100)


nPermut = 10000

permut(epi_trans, paste0("000-data/009-fdr_matched_libs/trans",job,".Rdata"))
permut(epi_cis, paste0("000-data/009-fdr_matched_libs/cis",job,".Rdata"))
