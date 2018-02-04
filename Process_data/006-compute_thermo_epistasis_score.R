library(parallel)
source("../000-functions.R")


#####################################################################################################################
### Functions  ######################################################################################################
#####################################################################################################################


# Function to compute FDR by permutation 
permut = function(d){
  
  set.seed(25061985)
  
  d = d[order(d$p),]
  nd = nrow(d)
  rp = d$p
  
  ep = as.matrix(d[,c("epi1","epi2","epi3")])
  
  # for each randomization
  # (need to split it in chunks of 100 because of memory issues)
  q = do.call("cbind",lapply(1:(nPermut/100),function(x){
    cat("\r", x, "\t\t\t")
    q = do.call("cbind",mclapply(1:100,function(y){
      
      # randomize
      ep = apply(ep,2,sample)
      
      # compute mean and standard error
      m = rowMeans(ep)
      s = sqrt(rowSums((ep-m)^2)/2) / sqrt(3)
      
      # t-test
      p = one.sample.ttest(m,s,0,3)
      
      # count number of false discoveries at each p-value threshold
      q = cut(p,breaks=c(0,rp))
      
      # cumulative sum to count all the false discoveries below a given p-value threshold
      cumsum(table(q))
      
    },mc.cores=8))
    
    m = rowMeans(q)
    s = sqrt(rowSums((q-m)^2)/(nPermut-1) / nPermut)
    
    cbind(m,s)
    
  }))
  
  # average number of false discoveries and standard error at each p-value threshold
  m = rowMeans(q[,seq(1,nPermut/100*2,2)])
  s = sqrt(rowSums(q[,seq(2,nPermut/100*2,2)]^2))
  
  # number of total discoveries at each p-value threshold
  o = 1:length(rp)
  
  # false discovery rate
  d$fdr = m/o
  # standard error of FDR is standard error of false discoveries divided by total number of discoveries
  d$fdrsd = s / o
  # because it is stochastic for the lowest p-values (No of false discoveries is low), it is better to define the p-value threshold of a given FDR as the maximal p-value for which FDR < 0.2
  # it is better to flag these hits instead of having to do this all the time
  d$hit = d$p <= max(rp[d$fdr < 0.2])
  
  d
}



# main function to compute epistasis scores
epi_scores = function(input_file, parameters, output_file){
  
  # load the data from the multiplicative epistasis score because signles and doubles are already matched
  # epistasis scores will be replaced by the ones calculated from the thermodynamic model
  load(input_file)
  
  # compute epistasis score for all replicates
  epi = do.call("cbind",mclapply(1:3,function(i){
    pr = thermo_model_pred(d[,paste0("s1_ppi",i)],d[,paste0("s2_ppi",i)],parameters[1],parameters[2],parameters[3])
    d[,paste0("do_ppi",i)] - pr
  },mc.cores=3))
  epi = cbind(epi,rowMeans(epi))
  d[,grep("epi",colnames(d))] = epi
  stderr = sqrt(rowSums((epi[,1:3]-epi[,4])^2)/2) / sqrt(3)
  d$p = one.sample.ttest(epi[,4],stderr,0,3)
  
  d = permut(d)
  
  save(d,file=output_file)
}


#####################################################################################################################
### Main  ###########################################################################################################
#####################################################################################################################


dir.create("000-data/006-epistasis_thermo_model", showWarnings = F)


# get thermodynamic model parameters
load("000-data/004-fitted_parameters.Rdata")


nPermut = 10000
epi_scores("000-data/003-mult_epi/trans.Rdata", p$trans, "000-data/006-epistasis_thermo_model/trans.Rdata")
print("lib1\n")
epi_scores("000-data/003-mult_epi/trans_restricted.Rdata", p$cis.trans, "000-data/006-epistasis_thermo_model/trans_restricted.Rdata")
print("lib2\n")
epi_scores("000-data/003-mult_epi/cis_restricted.Rdata", p$cis.trans, "000-data/006-epistasis_thermo_model/cis_restricted.Rdata")
print("lib3\n")
epi_scores("000-data/003-mult_epi/trans_filter100.Rdata", p$trans, "000-data/006-epistasis_thermo_model/trans_filter100.Rdata")
print("lib4\n")
epi_scores("000-data/003-mult_epi/trans_restricted_filter100.Rdata", p$cis.trans, "000-data/006-epistasis_thermo_model/trans_restricted_filter100.Rdata")
print("lib5\n")
epi_scores("000-data/003-mult_epi/cis_restricted_filter100.Rdata", p$cis.trans, "000-data/006-epistasis_thermo_model/cis_restricted_filter100.Rdata")



