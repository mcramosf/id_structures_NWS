
# Set Up --------------------------------------------------------

# loading packages
library(blockmodeling)
library(sna)
library(tidyverse)
library(doParallel)
library(doRNG)

# setting seed
set.seed(2024)

# Data Import -------------------------------------------------------------

# Identity Interference Network  

# Each cell represents the percentage of respondents holding
# the two identities for whom one identity (row) 
# interfered with the other (column).

interf_net <- readRDS("data/interf_net.RData")

# Blockmodeling Interfering Networks -------------------------------------

# 1. Pre-specified valued -------------------------------------------------


# Two blocks

# 1. core of burdensome identities (stigma and burdensome)
# 2. periphery of everything else

#     1    2
# 1   c    c
# 2   n    n

# 1. core of burdensome identities (stigma and burdensome)
# 2. periphery of everything else

#     1    2
# 1   c    c
# 2   n    n

# Three blocks

# 1. core of stigmatized and taxing identities
# 2. semicore of roles
# 3. periphery of everything else

#     1    2    3
# 1   c    c    c
# 2   c    n    n
# 3   n    n    n 

# Four blocks

# 1. core of stigmatized and taxing identities
# 2. semicore of roles
# 3. liberal
# 4. conservative

#     1    2    3   4
# 1   c    c    c   c
# 2   c    n    n   n
# 3   n    n    n   c 
# 4   n    n    c   n

# Five blocks
# This is the hypothesized full model

# 1. core of stigmatized identities
# 2. semicore of taxing identities
# 3. semicore of roles
# 4. liberal
# 5. conservative

#     1    2    3   4   5
# 1   c    c    c   c   c
# 2   c    c    c   n   n
# 3   c    c    n   n   n   
# 4   n    n    n   n   c 
# 5   n    n    n   c   n

# Six blocks or more
# Beyond the hypothesized full model

# 1. core of stigmatized identities
# 2. semicore of taxing identities
# 3. semicore of roles
# 4. liberal
# 5. conservative

#     1    2    3   4   5   6
# 1   c    c    c   c   c   c,n
# 2   c    c    c   n   n   c,n
# 3   c    c    n   n   n   c,n
# 4   n    n    n   n   c   c,n
# 5   n    n    n   c   n   c,n
# 6   c,n  c,n  c,n c,n c,n c,n

# Pre-specifying blockmodel and running for loop for different k

hyp_bm <- function(nbl){
  neg_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # first row is always complete
  neg_net_bm[1,1,1:nbl] <- "com"
  
  
  # second row
  # this is block dependent
  if(nbl == 2){
    neg_net_bm[1,2,1:2] <- "nul"
  } else if(nbl == 3){
    neg_net_bm[1,2,1] <- "com"
    neg_net_bm[1,2,2:3] <- "nul"
  } else if (nbl == 4){
    neg_net_bm[1,2,1] <- "com"
    neg_net_bm[1,2,2:4] <- "nul"
  } else if (nbl > 4){
    neg_net_bm[1,2,1:3] <- "com"
    neg_net_bm[1,2,4:5] <- "nul"
  }
  
  # third row
  if(nbl == 3){
    neg_net_bm[1,3,1:3] <- "nul"
  } else if (nbl == 4){
    neg_net_bm[1,3,1:3] <- "nul"
    neg_net_bm[1,3,4] <- "com"
  } else if (nbl > 4){
    neg_net_bm[1,3,1:2] <- "com"
    neg_net_bm[1,3,3:5] <- "nul"
  }
  
  # fourth row
  if (nbl == 4){
    neg_net_bm[1,4,1:2] <- "nul"
    neg_net_bm[1,4,3] <- "com"
    neg_net_bm[1,4,4] <- "nul"
  } else if (nbl > 4){
    neg_net_bm[1,4,1:4] <- "nul"
    neg_net_bm[1,4,5] <- "com"
  }
  
  #fifth row
  if (nbl > 4 ){
    neg_net_bm[1,5,1:3] <- "nul"
    neg_net_bm[1,5,4] <- "com"
    neg_net_bm[1,5,5] <- "nul"
  }
  
  # allowing cols to be complete or nul 
  #for 6th or greater cols
  if(nbl > 5){
    neg_net_bm[1,1:nbl,6:nbl] <- "com"
    neg_net_bm[2,1:nbl,6:nbl] <- "nul"
  }
  
  # allowing rows to be complete or nul 
  #for 6th or greater rows
  if(nbl > 5){
    neg_net_bm[1,6:nbl,1:5] <- "com"
    neg_net_bm[2,6:nbl,1:5] <- "nul"
  }
  
  return(neg_net_bm)
}

printIMmodel<-function(IM)print(apply(IM,2:3,function(x)paste(na.omit(x),collapse=",")))
printIMmodel(hyp_bm(5))

# Fitting the Specification
res_pres_neg_bm<-list()
res_pres_neg_dens<-list()

for(k in 2:10){
  res_pres_neg_bm[[as.character(k)]]<-optRandomParC(M=interf_net,
                                                    k=k,
                                                    preSpecM=.50,
                                                    rep=ifelse(k==1,1,1000),
                                                    blocks=hyp_bm(k),
                                                    approach="val",
                                                    nCores=ifelse(k==1,1,0))
  
  res_pres_neg_dens[[as.character(k)]]<-funByBlocks(interf_net,
                                                    clu=orderClu(res_pres_neg_bm[[as.character(k)]]),
                                                    na.rm=TRUE)
  res_pres_neg_dens[[as.character(k)]][is.nan(res_pres_neg_dens[[as.character(k)]])] <- 0
  
}

# Exploring the fit

res_pres_neg_bm[[5]]

#Saving these results so that we have something to work on later
saveRDS(res_pres_neg_bm,
        "results/res_pres_neg_bm.Rdata")

saveRDS(res_pres_neg_dens,
        "results/res_pres_neg_dens.Rdata")

# Relative Fit Function
rf_pre_neg <- list()

for (k in 1:9){
  rf_pre_neg[[k]] <- RF(res = res_pres_neg_bm[[k]], 
                        m = 30, 
                        loops = TRUE)
}

rf_pre_neg

# saveRDS(rf_pre_neg,
#         "manuscript/network_science/pre_proofs/reproducibility_files/results/res_pres_neg_rf.Rdata")

# 2. Cohesive  ---------------------------------------------

# Example with four blocks
#      1        2          3        4      
# 1   com      nul        nul      nul   
# 2   nul      com        nul      nul    
# 3   nul      nul        com      nul    
# 4   nul      nul        nul      com   

mcp <- function(nbl){
  neg_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # turning everything to null
  neg_net_bm[1,,1:nbl] <- "nul"
  
  for (i in 1:nbl){
    neg_net_bm[1,i,i] <- "com"
  }
  
  return(neg_net_bm)
}

# Fitting Model 

res_coh_neg_bm <- list()

for (k in 2:10){
  res_coh_neg_bm[[as.character(k)]] <- optRandomParC(M=interf_net,
                                                     k=k,
                                                     preSpecM=.50,
                                                     rep=ifelse(k==1,1,1000),
                                                     blocks=mcp(k),
                                                     approach="val",
                                                     nCores=ifelse(k==1,1,0))
}


#Saving results
saveRDS(res_coh_neg_bm,
        "results/res_coh_neg_bm.Rdata")

# Relative Fit Function
rf_coh_neg <- list()

for (k in 1:9){
  rf_coh_neg[[k]] <- RF(res = res_coh_neg_bm[[k]], 
                        m = 30, 
                        loops = TRUE)
}

rf_coh_neg

saveRDS(rf_coh_neg,
        "results/res_coh_neg_rf.Rdata")


# 3. Hierarchical ------------------------------------------------

#With a hierarchical blockmodel, 
#the nodes within clusters are not 
#linked to each other. The clusters in a
#hierarchical blockmodel can be ordered 
#sequentially in such a way that the nodes
#from each cluster (except for the highest-ranking clusters)
#are linked to the nodes from the cluster immediately

# Example with four blocks
#      1        2          3        4      
# 1   nul      nul        nul      nul   
# 2   com      nul        nul      nul    
# 3   nul      com        nul      nul    
# 4   nul      nul        com      nul   

mcp <- function(nbl){
  neg_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # turning everything to null
  neg_net_bm[1,,1:nbl] <- "nul"
  
  # 
  for (i in 1:(nbl-1)){
    neg_net_bm[1,(i+1),i] <- "com"
  }
  
  return(neg_net_bm)
}

printIMmodel(mcp(10))

# Fitting the Blockmodel
res_hier_neg_bm <- list()

for (k in 2:10){
  res_hier_neg_bm[[as.character(k)]] <- optRandomParC(M=interf_net,
                                                      k=k,
                                                      preSpecM=.50,
                                                      rep=ifelse(k==1,1,1000),
                                                      blocks=mcp(k),
                                                      approach="val",
                                                      nCores=ifelse(k==1,1,0))
}

#Saving these results

saveRDS(res_hier_neg_bm,
        "results/res_hier_neg_bm.Rdata")


# Relative Fit Function
rf_hier_neg <- list()

for (k in 1:9){
  rf_hier_neg[[k]] <- RF(res = res_hier_neg_bm[[k]], 
                         m = 30, 
                         loops = TRUE)
}

rf_hier_neg

saveRDS(rf_hier_neg,
        "results/res_hier_neg_rf.Rdata")

# 4. Hierarchical cohesive ------------------------------------------------
#same as hierarchical but now cohesive blocls

mcp <- function(nbl){
  neg_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # turning everything to null
  neg_net_bm[1,,1:nbl] <- "nul"
  
  # 
  for (i in 1:nbl){
    neg_net_bm[1,i,i] <- "com"
  }
  
  for (j in 1:(nbl-1)){
    neg_net_bm[1,(j+1),j] <- "com"
  }
  
  return(neg_net_bm)
}

printIMmodel(mcp(10))

# fitting the blockmodel

res_hier_coh_neg_bm <- list()

for (k in 2:10){
  res_hier_coh_neg_bm[[as.character(k)]] <- optRandomParC(M=interf_net,
                                                          k=k,
                                                          preSpecM=.50,
                                                          rep=ifelse(k==1,1,1000),
                                                          blocks=mcp(k),
                                                          approach="val",
                                                          nCores=ifelse(k==1,1,0))
}

#Saving these results 

saveRDS(res_hier_coh_neg_bm,
        "results/res_hier_coh_neg_bm.Rdata")

# Relative Fit Function
rf_hier_coh_neg <- list()

for (k in 1:9){
  rf_hier_coh_neg[[k]] <- RF(res = res_hier_coh_neg_bm[[k]], 
                             m = 30, 
                             loops = TRUE)
}

rf_hier_coh_neg

saveRDS(rf_hier_coh_neg,
        "results/res_hier_coh_neg_rf.Rdata")

# 5. Transitive  ------------------------------------------------

#A transitivity blockmodel is similar to a 
# hierarchical blockmodel. The exception is 
# when links exist from clusters on the lower
# level to all clusters on the higher levels. 
# A transitive-cohesive blockmodel is like a 
#transitivity blockmodel, but with the former 
#one links are found between nodes from the same cluster.

mcp <- function(nbl){
  neg_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # turning everything to null
  neg_net_bm[1,,1:nbl] <- "nul"
  
  # 
  for (i in 1:(nbl-1)){
    neg_net_bm[1,(i+1):nbl,i] <- "com"
  }
  
  return(neg_net_bm)
}

printIMmodel(mcp(10))

# fitting the blockmodel
res_trans_neg_bm <- list()

for (k in 2:10){
  res_trans_neg_bm[[as.character(k)]] <- optRandomParC(M=interf_net,
                                                       k=k,
                                                       preSpecM=.50,
                                                       rep=ifelse(k==1,1,1000),
                                                       blocks=mcp(k),
                                                       approach="val",
                                                       nCores=ifelse(k==1,1,0))
}

#Saving these results so that we have something to work on later
saveRDS(res_trans_neg_bm,
        "results/res_trans_neg_bm.Rdata")

# Relative Fit Function
rf_trans_neg <- list()

for (k in 1:9){
  rf_trans_neg[[k]] <- RF(res = res_trans_neg_bm[[k]], 
                          m = 30, 
                          loops = TRUE)
}

rf_trans_neg

saveRDS(rf_trans_neg,
        "results/res_trans_neg_rf.Rdata")


# 6. Transitive cohesive  ------------------------------------------------

#same as transitive but now cohesive blocks
mcp <- function(nbl){
  neg_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # turning everything to null
  neg_net_bm[1,,1:nbl] <- "nul"
  
  for (i in 1:nbl){
    neg_net_bm[1,i,i] <- "com"
  }
  
  for (i in 1:(nbl-1)){
    neg_net_bm[1,(i+1):nbl,i] <- "com"
  }
  
  return(neg_net_bm)
}

printIMmodel(mcp(10))

# fitting the blockmodel
res_trans_coh_neg_bm <- list()

for (k in 2:10){
  res_trans_coh_neg_bm[[as.character(k)]] <- optRandomParC(M=interf_net,
                                                           k=k,
                                                           preSpecM=.50,
                                                           rep=ifelse(k==1,1,1000),
                                                           blocks=mcp(k),
                                                           approach="val",
                                                           nCores=ifelse(k==1,1,0))
}

#Saving these results so that we have something to work on later
saveRDS(res_trans_coh_neg_bm,
        "results/res_trans_coh_neg_bm.Rdata")


# Relative Fit Function
rf_trans_coh_neg <- list()

for (k in 1:9){
  rf_trans_coh_neg[[k]] <- RF(res = res_trans_coh_neg_bm[[k]], 
                              m = 30, 
                              loops = TRUE)
}

rf_trans_coh_neg

saveRDS(rf_trans_coh_neg,
        "results/res_trans_coh_neg_rf.Rdata")

# The manuscript figures are obtained in the next script


