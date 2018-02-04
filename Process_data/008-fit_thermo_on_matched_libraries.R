job = as.numeric(as.character(commandArgs(trailingOnly = TRUE)))


#####################################################################################################################
### Load data and pool replicates ###################################################################################
#####################################################################################################################


load(paste0("000-data/007-matched_libs/trans",job,".Rdata"))
load(paste0("000-data/007-matched_libs/cis",job,".Rdata"))


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

  


#####################################################################################################################
### Fit model  ######################################################################################################
#####################################################################################################################

# get the thermodynamic model function that predicts double mutant PPI scores from single mutants
source("../000-functions.R")


# parameters to search
sx = seq(1.1,1.6,0.01)
sy = seq(1.1,1.6,0.01)
sb = seq(0.3,1,0.01)



# fit model on trans
rss = do.call("c",lapply(as.list(sx),function(x){
  do.call("c",lapply(as.list(sy),function(y){
    if(y < x){
        rep(Inf,length(sb))
    }else{
      do.call("c",lapply(as.list(sb),function(B){
        # prediction from the model
        pr = thermo_model_pred(tr$s1,tr$s2,x,y,B)
        # residual sum of square
        sum((tr$do-pr)^2)
      }))
    }
  }))
}))

rss = array(rss, dim=c(length(sb),length(sy),length(sx)))

min.rss.tr = min(rss)

# retrieve the best parameters
w = which(rss == min.rss.tr, arr.ind=T)
x.trans = sx[w[1,3]]
y.trans = sy[w[1,2]]
B.trans = sb[w[1,1]]



# fit model on cis
rss = do.call("c",lapply(as.list(sx),function(x){
    do.call("c",lapply(as.list(sy),function(y){
        if(y < x){
            rep(Inf,length(sb))
        }else{
            do.call("c",lapply(as.list(sb),function(B){
                # prediction from the model
                pr = thermo_model_pred(ci$s1,ci$s2,x,y,B)
                # residual sum of square
                sum((ci$do-pr)^2)
            }))
        }
    }))
}))

rss = array(rss, dim=c(length(sb),length(sy),length(sx)))

min.rss.ci = min(rss)

# retrieve the best parameters
w = which(rss == min.rss.ci, arr.ind=T)
x.cis = sx[w[1,3]]
y.cis = sy[w[1,2]]
B.cis = sb[w[1,1]]



# evaluate model fitted on trans on the cis library
pr = thermo_model_pred(ci$s1,ci$s2,x.trans,y.trans,B.trans)
rss.tr.ci = sum((ci$do-pr)^2)
# evaluate model fitted on cis on the trans library
pr = thermo_model_pred(tr$s1,tr$s2,x.cis,y.cis,B.cis)
rss.ci.tr = sum((tr$do-pr)^2)



var_exp.tr.tr = 1 - min.rss.tr/sum((tr$do - mean(tr$do))^2)
var_exp.tr.ci = 1 - rss.tr.ci/sum((ci$do - mean(ci$do))^2)
var_exp.ci.ci = 1 - min.rss.ci/sum((ci$do - mean(ci$do))^2)
var_exp.ci.tr = 1 - rss.ci.tr/sum((tr$do - mean(tr$do))^2)

rmsd.tr.tr = sqrt(min.rss.tr/nrow(tr))
rmsd.tr.ci = sqrt(rss.tr.ci/nrow(ci))
rmsd.ci.ci = sqrt(min.rss.ci/nrow(ci))
rmsd.ci.tr = sqrt(rss.ci.tr/nrow(tr))


out = c(x.trans,y.trans,B.trans,var_exp.tr.tr,var_exp.tr.ci,rmsd.tr.tr,rmsd.tr.ci,   # model fitted on trans
        x.cis,y.cis,B.cis,var_exp.ci.ci,var_exp.ci.tr,rmsd.ci.ci,rmsd.ci.tr          # model fitted on cis
)


save(out,file=paste0("000-data/008-fit_thermo_on_matched_libraries/",job,".Rdata"))

