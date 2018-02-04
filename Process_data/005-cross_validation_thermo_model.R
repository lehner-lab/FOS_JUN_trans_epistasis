
job = as.numeric(as.character(commandArgs(trailingOnly = TRUE)))
lib = ceiling(job/1000)

set.seed(job)

#####################################################################################################################
### Load data, pool replicates and randomly split the data in two independent quarters based on singles mutants id  #
#####################################################################################################################

# load data from the multiplicative model because double mutants are already matched to singles
# epistasis arenot going to be used here anyway
input_files = c(
  "000-data/003-mult_epi/trans.Rdata",
  "000-data/003-mult_epi/trans_restricted.Rdata",
  "000-data/003-mult_epi/cis_restricted.Rdata"
)


# if it is the trans library
if(lib <= 2){
  
  load(input_files[lib])

  # randomize single mutants id
  rand_id1 = sample(unique(d$id1))
  rand_id2 = sample(unique(d$id2))
  
  # split single mutants into test and train sets
  id1_test = rand_id1[1:(length(rand_id1)/2)]
  id2_test = rand_id2[1:(length(rand_id2)/2)]
  id1_train = rand_id1[!(rand_id1 %in% id1_test)]
  id2_train = rand_id2[!(rand_id2 %in% id2_test)]
  
  # pool replicates
  d = data.frame(
    id1 = rep(d$id1,3),
    id2 = rep(d$id2,3),
    s1 = c(d$s1_ppi1,d$s1_ppi2,d$s1_ppi3),
    s2 = c(d$s2_ppi1,d$s2_ppi2,d$s2_ppi3),
    do = c(d$do_ppi1,d$do_ppi2,d$do_ppi3)
  )
  
  # split double mutants into train and test sets
  train = d[d$id1 %in% id1_train & d$id2 %in% id2_train,]
  test = d[d$id1 %in% id1_test & d$id2 %in% id2_test,]
    
}else if(lib == 3){
    
    load(input_files[lib])
    
    # randomize single mutants id
    rand_id1 = sample(unique(c(d$id1,d$id2)))
    id1_test = rand_id1[1:(length(rand_id1)/2)]
    id1_train = rand_id1[!(rand_id1 %in% id1_test)]

    # pool replicates
    d = data.frame(
        id1 = rep(d$id1,3),
        id2 = rep(d$id2,3),
        s1 = c(d$s1_ppi1,d$s1_ppi2,d$s1_ppi3),
        s2 = c(d$s2_ppi1,d$s2_ppi2,d$s2_ppi3),
        do = c(d$do_ppi1,d$do_ppi2,d$do_ppi3)
    )

    # split double mutants into train and test sets
    train = d[d$id1 %in% id1_train & d$id2 %in% id1_train,]
    test = d[d$id1 %in% id1_test & d$id2 %in% id1_test,]
    
    
}else{  # if lib == 4, the model is fitted on the cis and trans restricted libs pooled together
  
  load(input_files[2])
  tr = d
  load(input_files[3])
  ci = d
  
  # randomize single mutants id and split them into test and train sets
  # Fos single mutants will be split in the same way between the two libraries
  rand_id1 = sample(unique(c(ci$id1,ci$id2,tr$id1)))
  rand_id2 = sample(unique(tr$id2))
  id1_test = rand_id1[1:(length(rand_id1)/2)]
  id2_test = rand_id2[1:(length(rand_id2)/2)]
  id1_train = rand_id1[!(rand_id1 %in% id1_test)]
  id2_train = rand_id2[!(rand_id2 %in% id2_test)]
  
  # pool replicates
  ci = data.frame(
    id1 = rep(ci$id1,3),
    id2 = rep(ci$id2,3),
    s1 = c(ci$s1_ppi1,ci$s1_ppi2,ci$s1_ppi3),
    s2 = c(ci$s2_ppi1,ci$s2_ppi2,ci$s2_ppi3),
    do = c(ci$do_ppi1,ci$do_ppi2,ci$do_ppi3)
  )
  tr = data.frame(
    id1 = rep(tr$id1,3),
    id2 = rep(tr$id2,3),
    s1 = c(tr$s1_ppi1,tr$s1_ppi2,tr$s1_ppi3),
    s2 = c(tr$s2_ppi1,tr$s2_ppi2,tr$s2_ppi3),
    do = c(tr$do_ppi1,tr$do_ppi2,tr$do_ppi3)
  )

  # split double mutants into train and test sets
  train.ci = ci[ci$id1 %in% id1_train & ci$id2 %in% id1_train,]
  test.ci = ci[ci$id1 %in% id1_test & ci$id2 %in% id1_test,]
  train.tr = tr[tr$id1 %in% id1_train & tr$id2 %in% id2_train,]
  test.tr = tr[tr$id1 %in% id1_test & tr$id2 %in% id2_test,]

  # pool the two libraries
  train = rbind(train.ci, train.tr)  
  test = rbind(test.ci, test.tr)  
}
  


#####################################################################################################################
### Fit model  ######################################################################################################
#####################################################################################################################


# get the thermodynamic model function that predicts double mutant PPI scores from single mutants
source("../000-functions.R")


# parameters to search
sx = seq(1.1,1.6,0.01)
sy = seq(1.1,1.6,0.01)
sb = seq(0.3,1,0.01)

# fit model on train set
rss = do.call("c",lapply(as.list(sx),function(x){
  do.call("c",lapply(as.list(sy),function(y){
    if(y < x){
      rep(Inf,length(sb))
    }else{
      do.call("c",lapply(as.list(sb),function(B){
        # prediction from the model
        pr = thermo_model_pred(train$s1,train$s2,x,y,B)
        # residual sum of square
        sum((train$do-pr)^2)
      }))
    }
  }))
}))

rss = array(rss, dim=c(length(sb),length(sy),length(sx)))


# retrieve the best parameters
w = which(rss == min(rss), arr.ind=T)
x = sx[w[1,3]]
y = sy[w[1,2]]
B = sb[w[1,1]]



# assess fitted model on test set
pr = thermo_model_pred(test$s1,test$s2,x,y,B)
rss_test = sum((test$do-pr)^2)

var_exp_train = 1 - min(rss)/sum((train$do - mean(train$do))^2)
var_exp_test = 1 - rss_test/sum((test$do - mean(test$do))^2)
rmsd = sqrt(rss_test/nrow(test))

out = c(x,y,B,var_exp_train,var_exp_test,rmsd)


save(out,file=paste0("000-data/005-CV_thermo/",job,".Rdata"))