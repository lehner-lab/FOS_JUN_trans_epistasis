library(parallel)
library(stringdist)



#####################################################################################################################
### Load and format raw data ########################################################################################
#####################################################################################################################

prep_data_trans = function(input_file){
    
    d = read.delim(input_file)
    
    # convert NAs for wt and single variants to -1 (pos) or wt (aa)
    d$aa_pos_F[d$F %in% c("wt","S")] = -1
    d$aa_pos_J[d$J %in% c("wt","S")] = -1
    d$wt_aa_F[d$F %in% c("wt","S")] = "wt"
    d$wt_aa_J[d$J %in% c("wt","S")] = "wt"
    d$mut_aa_F[d$F %in% c("wt","S")] = "wt"
    d$mut_aa_J[d$J %in% c("wt","S")] = "wt"
    
    
    # aggregate to aa
    d = aggregate(d[,4:9],as.list(d[,c(13:15,19:21)]),sum)
    names(d)[1:6] = c("pos1","wt1","mut1","pos2","wt2","mut2")
    d = d[,c(1,4,2,5,3,6,7:12)]
    
    # convert position from 0:6 to 1:7
    d$pos1 = d$pos1 + 1
    d$pos2 = d$pos2 + 1
    
    # variant IDs
    d = cbind(
        id1 = paste(ceiling(d$pos1/7),letters[(d$pos1-1)%%7+1],d$mut1,sep=""),
        id2 = paste(ceiling(d$pos2/7),letters[(d$pos2-1)%%7+1],d$mut2,sep=""),
        d
    )
}

prep_data_cis = function(input_file){
    
    d = read.delim(input_file)
    
    # total read count before filtering out variants with more than 2 aa substitutions
    tot_read_count = colSums(d[,3:8])
    
    # filter out variants with more than 2 aa substitutions
    pos = strsplit(d$ns_aa_pos, " ")
    nb = do.call("c",mclapply(pos,length,mc.cores=4))
    d = d[nb <= 2,]
    
    # convert NAs for wt and single variants to 0 (pos) or wt (aa)
    d$ns_aa_pos[d$F %in% c("wt","S")] = -1
    d$ns_aa_wt[d$F %in% c("wt","S")] = "wt"
    d$ns_aa_mut[d$F %in% c("wt","S")] = "wt"
    
    
    # aggregate to aa
    d = aggregate(d[,3:8],as.list(d[,15:17]),sum)
    
    # format the dataframe
    pos = strsplit(d$ns_aa_pos, " ")
    nb = do.call("c",mclapply(pos,length,mc.cores=4))
    doubles = d[nb == 2,]
    pos = do.call("rbind",pos[nb == 2])
    wt = do.call("rbind",strsplit(doubles$ns_aa_wt, " "))
    mut = do.call("rbind",strsplit(doubles$ns_aa_mut, " "))
    doubles = cbind(pos,wt,mut,doubles[,4:9])
    names(doubles)[1:6] = c("pos1","pos2","wt1","wt2","mut1","mut2")
    
    singles = d[nb < 2,] # including wt
    singles$pos1 = -1
    singles$wt1 = "wt"
    singles$mut1 = "wt"
    singles = singles[,c(10,1,11,2,12,3,4:9)]
    names(singles)[1:6] = names(doubles)[1:6]
    
    d = rbind(singles,doubles)
    d$pos1 = as.numeric(d$pos1) + 1
    d$pos2 = as.numeric(d$pos2) + 1
    
    
    # variant IDs
    d = cbind(
        id1 = paste(ceiling(d$pos1/7),letters[(d$pos1-1)%%7+1],d$mut1,sep=""),
        id2 = paste(ceiling(d$pos2/7),letters[(d$pos2-1)%%7+1],d$mut2,sep=""),
        d
    )
    
    list(d, tot_read_count)
}


input_files_trans = list(
    "000-data/001-raw_data/trans_Q20.txt",
    "000-data/001-raw_data/trans_Q30.txt",
    "000-data/001-raw_data/trans_Q35.txt",
    "000-data/001-raw_data/trans_mutQualFiltered_Q20.txt"
)


input_files_cis = list(
    "000-data/001-raw_data/cis_Q20.txt",
    "000-data/001-raw_data/cis_Q30.txt",
    "000-data/001-raw_data/cis_Q35.txt",
    "000-data/001-raw_data/cis_mutQualFiltered_Q20.txt"
)



trans = mclapply(input_files_trans,function(x){
    prep_data_trans(x)
},mc.cores=4)

cis = mclapply(input_files_cis,function(x){
    prep_data_cis(x)
},mc.cores=4)





#####################################################################################################################
### Compute PPI scores ##############################################################################################
#####################################################################################################################

# Function to compute PPI scores

ppi_scores = function(d, type_lib, OD, out_file){
    
    if(type_lib == "trans"){
        freq = t(d[,9:14])/colSums(d[,9:14])
    }else{
        tot_read = d[[2]]
        d = d[[1]]
        freq = t(d[,9:14])/tot_read
    }
    gen = log2(freq[4:6,] / freq[1:3,] * OD)
    wt = gen[,d$pos1 == 0 & d$pos2 == 0]
    
    ppi = t(gen / wt)
    mppi = rowMeans(ppi)
    
    d = cbind(d,ppi,mppi)
    names(d)[15:17] = c("ppi1","ppi2","ppi3")
    
    save(d,file=out_file)
    
    d
}



# ODs during competition assay
OD.ratio.trans = c(1.57,1.6,1.53) / rep(0.05,3) * c(2.66,2.49,2.55)  / rep(0.05,3)
OD.ratio.cis = c(1.585,1.615,1.65) / rep(0.1,3) * c(1.475,1.59,1.49)  / rep(0.05,3)


# output files
dir.create("000-data/002-PPIscores", showWarnings = F, recursive = T)
out_files = c(
    "000-data/002-PPIscores/trans_Q20.Rdata",
    "000-data/002-PPIscores/trans_Q30.Rdata",
    "000-data/002-PPIscores/trans_Q35.Rdata",
    "000-data/002-PPIscores/trans_mutQualFiltered_Q20.Rdata",
    "000-data/002-PPIscores/cis_Q20.Rdata",
    "000-data/002-PPIscores/cis_Q30.Rdata",
    "000-data/002-PPIscores/cis_Q35.Rdata",
    "000-data/002-PPIscores/cis_mutQualFiltered_Q20.Rdata"
)

# main
l = do.call("c",mclapply(1:4,function(i){
    list(   ppi_scores(trans[[i]], "trans", OD.ratio.trans, out_files[i]),
            ppi_scores(cis[[i]], "cis", OD.ratio.cis, out_files[i+4])
    )
}, mc.cores=4))



#####################################################################################################################
### Filter based on read counts and align modes between replicates ##################################################
#####################################################################################################################

# function to calculate the mode of detrimental variants
mode = function(x1){
    x1 = density(x1)
    x1 = cbind(x1$x,x1$y)
    x1 = cbind(x1,NA)
    for(i in 2:nrow(x1)){
        x1[i,3] = x1[i,2] - x1[i-1,2]
    }
    out = NA
    for(i in 2:nrow(x1)){
        if(x1[i,3] < 0){
            break
        }else{
            out = x1[i,1]
        }
    }
    out
}


# function to filter and align modes
align = function(d, out_file){
    
    # filter based on read counts
    d = d[d$i1 > 10 & d$i2 > 10 & d$i3 > 10 & d$o1 > 0 & d$o2 > 0 & d$o3 > 0,]
    # filter out stops and the WT
    d = d[d$mut1 != "*" & d$mut2 != "*" & (d$pos1 != 0 | d$pos2 != 0),]
    
    # identify mode of detrimental variants
    m = apply(d[,c("ppi1","ppi2","ppi3")],2,mode)
    # correction mode
    m = (m - max(m))/(1 - max(m))
    # correction PPI scores
    d$ppi1 = (d$ppi1 - m[1]) / (1 - m[1])
    d$ppi2 = (d$ppi2 - m[2]) / (1 - m[2])
    d$ppi3 = (d$ppi3 - m[3]) / (1 - m[3])
    d$mppi = rowMeans(d[,c("ppi1","ppi2","ppi3")])
    
    save(d,file=out_file)
    
    d
}




# output files
out_files = c(
    "000-data/002-PPIscores/trans_Q20_readCountFiltered.Rdata",
    "000-data/002-PPIscores/trans_Q30_readCountFiltered.Rdata",
    "000-data/002-PPIscores/trans_Q35_readCountFiltered.Rdata",
    "000-data/002-PPIscores/trans_mutQualFiltered_Q20_readCountFiltered.Rdata",
    "000-data/002-PPIscores/cis_Q20_readCountFiltered.Rdata",
    "000-data/002-PPIscores/cis_Q30_readCountFiltered.Rdata",
    "000-data/002-PPIscores/cis_Q35_readCountFiltered.Rdata",
    "000-data/002-PPIscores/cis_mutQualFiltered_Q20_readCountFiltered.Rdata"
)

# sort the list of data frames to have them in the same order as output files
l = l[c(1,3,5,7,2,4,6,8)]


# main
l = mclapply(1:8,function(i){
    align(l[[i]],out_files[i])
},mc.cores=4)




#####################################################################################################################
### Filter out amino acid substitutions not reachable by a single aa change #########################################
### (for the cis/trans comparison in the second part of the paper)          #########################################
### Needs to be done for cis single mutants as well                         #########################################
#####################################################################################################################

# done only on Q20 without mutQualFiltered, because the other ones are used only to show that more stringent filtering minimally affects PPI score


# Minimal hamming distance between wt and mut codons
gc = read.delim("000-data/000-genetic_code.txt",header=F)
# add a dummy codon to deal with single mutants. Will be corrected afterward (cf below)
gc = rbind(gc,c("XXX","wt"))


restrict = function(d, type, out_file){
    # get the WT codon
    s1 = d$pos1 * 3 - 2
    s2 = d$pos2 * 3 - 2
    e1 = d$pos1 * 3
    e2 = d$pos2 * 3
    
    wtSeq1 = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA"
    if(type == "trans"){
        wtSeq2 = "ATCGCCCGGCTGGAGGAAAAAGTGAAAACCTTGAAAGCTCAGAACTCGGAGCTGGCGTCCACGGCCAACATGCTCAGGGAACAGGTGGCACAGCTT"
    }else{
        wtSeq2 = "ACTGATACACTCCAAGCGGAGACAGACCAACTAGAAGATGAGAAGTCTGCTTTGCAGACCGAGATTGCCAACCTGCTGAAGGAGAAGGAAAAACTA"
    }
    
    wt_cod1 = sapply(1:length(s1),function(x){
        substr(wtSeq1,s1[x],e1[x])
    })
    wt_cod2 = sapply(1:length(s2),function(x){
        substr(wtSeq2,s2[x],e2[x])
    })
    
    # calculate hamming distance between WT codon and closest codon encoding the mutant amino acid
    hamm = do.call("rbind",mclapply(1:nrow(d),function(x){
        mut_cod1 = gc$V1[gc$V2 == d$mut1[x]]
        mut_cod2 = gc$V1[gc$V2 == d$mut2[x]]
        hamm1 = min(stringdist(wt_cod1[x],mut_cod1))
        hamm2 = min(stringdist(wt_cod2[x],mut_cod2))
        c(hamm1,hamm2)
    },mc.cores=4))
    # this will have given a wrong result for single mutant, so it need to be corrected (more efficient like than this than an if statement in the sapply)
    hamm[,1][d$pos1 == 0] = 0
    hamm[,2][d$pos2 == 0] = 0
    
    
    d = d[hamm[,1] <= 1 & hamm[,2] <= 1,]
    
    save(d, file = out_file)
}


restrict(l[[1]], "trans", "000-data/002-PPIscores/trans_Q20_restricted.Rdata")
restrict(l[[5]], "cis","000-data/002-PPIscores/cis_Q20_restricted.Rdata")