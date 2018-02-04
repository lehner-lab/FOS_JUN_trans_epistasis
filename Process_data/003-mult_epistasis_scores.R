library(parallel)
source("../000-functions.R")



#####################################################################################################################
### Compute epistasis scores ########################################################################################
#####################################################################################################################


multEpi = function(d, type){
    
    # split data into sungles and doubles
    do = d[d$pos1 != 0 & d$pos2 != 0,]
    if(type == "trans"){
        s1 = d[d$pos2 == 0,c(1,3,5,7,9:18)]
        rownames(s1) = s1$id1
    }else{
        s1 = d[d$pos1 == 0,c(2,4,6,8:18)]
        rownames(s1) = s1$id2
        names(s1)[1:4] = c("id1","pos1","wt1","mut1")
    }
    s2 = d[d$pos1 == 0,c(2,4,6,8:18)]
    rownames(s2) = s2$id2

    # remove double mutants for which the single mutant is missing
    do = do[do$id1 %in% s1$id1 & do$id2 %in% s2$id2,]
    
    # match double mutants to the corresponding singles
    do = cbind(s1[do$id1,],s2[do$id2,],do[,9:18])
    
    # reorganize the data frame
    do = do[,c(1,15,2,16,3,17,4,18,5:10,19:24,29:34,11:14,25:28,35:38)]
    names(do)[c(21:26,35:38)] = paste0("do_",names(do)[c(9:14,27:30)])
    names(do)[c(15:20,31:34)] = paste0("s2_",names(do)[c(9:14,27:30)])
    names(do)[c(9:14,27:30)] = paste0("s1_",names(do)[c(9:14,27:30)])
    
    # compute epistasis scores
    do$epi1 = do$do_ppi1 - do$s1_ppi1 * do$s2_ppi1
    do$epi2 = do$do_ppi2 - do$s1_ppi2 * do$s2_ppi2
    do$epi3 = do$do_ppi3 - do$s1_ppi3 * do$s2_ppi3
    tmpMat = as.matrix(do[,c("epi1","epi2","epi3")])
    do$mepi = rowMeans(tmpMat)
    
    # test signficance
    stderr = sqrt(rowSums((tmpMat-do$mepi)^2)/2) / sqrt(3)
    do$p = one.sample.ttest(do$mepi,stderr,0,3)
    
    do
}


#####################################################################################################################
### Compute FDR by permutations #####################################################################################
#####################################################################################################################

permut = function(d, out_file){
  
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
  
  save(d,file=out_file)
}




#####################################################################################################################
### Main ############################################################################################################
#####################################################################################################################



input_files = c(
  "000-data/002-PPIscores/trans_Q20_readCountFiltered.Rdata",
  "000-data/002-PPIscores/trans_Q20_restricted.Rdata",             # restricted are already read count filtered
  "000-data/002-PPIscores/cis_Q20_restricted.Rdata"
)


dir.create("000-data/003-mult_epi", showWarnings = F, recursive = T)
output_files = c(
  "000-data/003-mult_epi/trans.Rdata",
  "000-data/003-mult_epi/trans_restricted.Rdata",
  "000-data/003-mult_epi/cis_restricted.Rdata"
)
output_files2 = c(
  "000-data/003-mult_epi/trans_filter100.Rdata",
  "000-data/003-mult_epi/trans_restricted_filter100.Rdata",
  "000-data/003-mult_epi/cis_restricted_filter100.Rdata"
)

type = c(
  "trans",
  "trans",
  "cis"
)

nPermut = 10000


sapply(1:3,function(lib){

    load(input_files[lib])
    d = multEpi(d, type[lib])
    permut(d, output_files[lib])
    
    d = d[d$s1_i1 > 100 & d$s1_i2 > 100 & d$s1_i3 > 100 & d$s2_i1 > 100 & d$s2_i2 > 100 & d$s2_i3 > 100 & d$do_i1 > 100 & d$do_i2 > 100 & d$do_i3 > 100,]
    permut(d, output_files2[lib])

})
