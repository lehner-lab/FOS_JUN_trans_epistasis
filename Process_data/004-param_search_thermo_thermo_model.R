
job = as.numeric(as.character(commandArgs(trailingOnly = TRUE)))
lib = ceiling(job/901)

#####################################################################################################################
### Load data and pool replicates  ##################################################################################
#####################################################################################################################

# load data from the multiplicative model because double mutants are already matched to singles
# epistasis arenot going to be used here anyway
input_files = c(
  "000-data/003-mult_epi/trans.Rdata",
  "000-data/003-mult_epi/trans_restricted.Rdata",
  "000-data/003-mult_epi/cis_restricted.Rdata"
)


# if lib == 4, the model is fitted on the cis and trans restricted libs pooled together
if(lib <= 3){

  load(input_files[lib])

  # pool replicates
  d = data.frame(
    s1 = c(d$s1_ppi1,d$s1_ppi2,d$s1_ppi3),
    s2 = c(d$s2_ppi1,d$s2_ppi2,d$s2_ppi3),
    do = c(d$do_ppi1,d$do_ppi2,d$do_ppi3)
  )

}else{
  
  load(input_files[2])
  tr = d
  load(input_files[3])
  ci = d
  
  # pool replicates
  ci = data.frame(
    s1 = c(ci$s1_ppi1,ci$s1_ppi2,ci$s1_ppi3),
    s2 = c(ci$s2_ppi1,ci$s2_ppi2,ci$s2_ppi3),
    do = c(ci$do_ppi1,ci$do_ppi2,ci$do_ppi3)
  )
  tr = data.frame(
    s1 = c(tr$s1_ppi1,tr$s1_ppi2,tr$s1_ppi3),
    s2 = c(tr$s2_ppi1,tr$s2_ppi2,tr$s2_ppi3),
    do = c(tr$do_ppi1,tr$do_ppi2,tr$do_ppi3)
  )

  d = rbind(ci,tr)
}
  



#####################################################################################################################
### Parameter search  ###############################################################################################
#####################################################################################################################

# get the thermodynamic model function that predicts double mutant PPI scores from single mutants
source("../000-functions.R")


# B parameter for this job
sb = seq(0.3,1.2,0.001)
B = sb[((job-1) %% length(sb)) + 1]

# Fit x and y for this value of B
sx = seq(1.1,2,0.01)
sy = seq(1.1,2,0.01)
rss = do.call("c",lapply(1:length(sx),function(i){
  x = sx[i]
  # the two parameters are equally intervertible, so no need to recalculate when y < x
  do.call("c",lapply(i:length(sy),function(j){
    y = sy[j]
    # prediction from the model
    pr = thermo_model_pred(d$s1,d$s2,x,y,B)
    # residual sum of square
    sum((d$do-pr)^2)
  }))
}))

save(rss,file=paste0("000-data/004-param_search/",job,".Rdata"))




