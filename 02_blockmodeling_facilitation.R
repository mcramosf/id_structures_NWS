
# Set Up --------------------------------------------------------

# loading packages
library(blockmodeling)
library(sna)
library(tidyverse)
library(doParallel)
library(doRNG)

# setting seed
set.seed(1234)

# Data Import -------------------------------------------------------------

# Identity Facilitation Network 

# Each cell represents the percentage of respondents holding
# the two identities for whom one identity (row) 
# facilitated the other (column).

facil_net <- readRDS("data/facil_net.RData")

# Identity Interference Network  

# Each cell represents the percentage of respondents holding
# the two identities for whom one identity (row) 
# interfered with the other (column).

interf_net <- readRDS("data/interf_net.RData")


# Blockmodeling Facilitation Networks -------------------------------------

# 1. Pre-specified valued -------------------------------------------------

# Pre-specified blockmodeling core-periphery and
# cohesive mix

# The hypothesized blockmodel has 5 positions. For block-
# models with fewer than 5 positions, I retain the core-
# periphery structure

# Two blocks 
# symmetric core

#     1    2
# 1   c    c
# 2   c    n

# Three blocks
# symmetric core 
#     1    2    3
# 1   c    c    c
# 2   c    c    c 
# 3   c    n    n 

# Four blocks
# includes periphery
#     1    2    3   4
# 1   c    c    c   c
# 2   c    c    c   n
# 3   c    c,n  c  c,n 
# 4   c    n    n   n

# Hypothesized Blockmodel
#Specification in Figure 1 (Panel A) of main manuscript
#      1        2          3        4     5 
# 1   com      com        com      com   c,n
# 2   com      com        com      com   c,n 
# 3   com      c,n        com      c,n   c,n    
# 4   com      nul        nul      nul   c,n   
# 5   c,n      c,n        c,n      c,n   c,n

# Setting Specification

hyp_bm <- function(nbl){
  pos_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # first row is complete for first 4 positions
  pos_net_bm[1,1,1:nbl] <- "com"
  if(nbl > 4){
    pos_net_bm[2,1,5:nbl] <- "nul"
  }
  
  # second row
  # this is block dependent
  if(nbl == 2){
    pos_net_bm[1,2,1] <- "com"
    pos_net_bm[1,2,2] <- "nul"
  } else if(nbl == 3){
    pos_net_bm[1,2,1:2] <- "com"
    pos_net_bm[1,2,3] <- "nul"
  } else if (nbl == 4){
    pos_net_bm[1,2,1:3] <- "com"
    pos_net_bm[1,2,4] <- "nul"
  } else if (nbl > 4){
    pos_net_bm[1,2,1:nbl] <- "com"
    pos_net_bm[2,2,5:nbl] <- "nul"
  }
  
  # third row
  if(nbl == 3){
    pos_net_bm[1,3,1] <- "com"
    pos_net_bm[1,3,2:3] <- "nul"
  } else if (nbl == 4){
    pos_net_bm[1,3,1:4] <- "com"
    pos_net_bm[2,3,c(2,4)] <- "nul"
  } else if (nbl > 4){
    pos_net_bm[1,3,1:nbl] <- "com"
    pos_net_bm[2,3,c(2,4:nbl)] <- "nul"
  }
  
  # fourth row
  if (nbl == 4){
    pos_net_bm[1,4,1] <- "com"
    pos_net_bm[1,4,2:4] <- "nul"
  } else if (nbl > 4){
    pos_net_bm[1,4,c(1,5:nbl)] <- "com"
    pos_net_bm[1,4,2:4] <- "nul"
    pos_net_bm[2,4,5:nbl] <- "nul"
  }
  
  #fifth row
  if (nbl > 4 ){
    pos_net_bm[1,5,1:nbl] <- "com"
    pos_net_bm[2,5,1:nbl] <- "nul"
  }
  
  # allowing cols to be complete or nul 
  #for 6th or greater cols
  if(nbl > 5){
    pos_net_bm[1,1:nbl,6:nbl] <- "com"
    pos_net_bm[2,1:nbl,6:nbl] <- "nul"
  }
  
  # allowing rows to be complete or nul 
  #for 6th or greater rows
  if(nbl > 5){
    pos_net_bm[1,6:nbl,1:5] <- "com"
    pos_net_bm[2,6:nbl,1:5] <- "nul"
  }
  
  return(pos_net_bm)
}

printIMmodel<-function(IM)print(apply(IM,2:3,function(x)paste(na.omit(x),collapse=",")))
# here you can check how the blockmodels look 
printIMmodel(hyp_bm(4))

# Fitting the Specification
res_pres_pos_bm<-list()
res_pres_pos_dens<-list()

for(k in 2:10){
  res_pres_pos_bm[[as.character(k)]]<-optRandomParC(M=facil_net,
                                                    k=k,
                                                    preSpecM=.50,
                                                    rep=ifelse(k==1,1,1000),
                                                    blocks=hyp_bm(k),
                                                    approach="val",
                                                    nCores=ifelse(k==1,1,0))
  
  res_pres_pos_dens[[as.character(k)]]<-funByBlocks(facil_net,
                                                    clu=orderClu(res_pres_pos_bm[[as.character(k)]]),
                                                    na.rm=TRUE)
  res_pres_pos_dens[[as.character(k)]][is.nan(res_pres_pos_dens[[as.character(k)]])] <- 0
  

}

# Saving these results for later use

# Blockmodeling 
saveRDS(res_pres_pos_bm,
        "results/res_pres_pos_bm.Rdata")

# Density Matrices
saveRDS(res_pres_pos_dens,
        "results/res_pres_pos_dens.Rdata")

# Exploring results
res_pres_pos_bm[[5]]
res_pres_pos_dens[[5]]

# Obtaining Relative Fit 
rf_pres_pos <- list()

for (k in 1:9){
  rf_pres_pos[[k]] <- RF(res = res_pres_pos_bm[[k]], 
                         m = 30, 
                         loops = TRUE)
}

rf_pres_pos

# saving results
saveRDS(rf_pres_pos,
        "results/res_pres_pos_rf.Rdata")

# let's take a look
map(rf_pres_pos, 
    ~ unlist(.x$RF))


# 2. Transitive -------------------------------------------------

#A transitivity blockmodel is similar to a 
# hierarchical blockmodel. The exception is 
# when links exist from clusters on the lower
# level to all clusters on the higher levels. 
# A transitive-cohesive blockmodel is like a 
# transitivity blockmodel, but with the former 
# one links are found between nodes from the same cluster.

mcp <- function(nbl){
  pos_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # turning everything to null
  pos_net_bm[1,,1:nbl] <- "nul"
  
  # 
  for (i in 1:(nbl-1)){
    pos_net_bm[1,(i+1):nbl,i] <- "com"
  }
  
  return(pos_net_bm)
}

printIMmodel(mcp(10))

# fitting the blockmodel

res_trans_pos_bm <- list()

for (k in 2:10){
  res_trans_pos_bm[[as.character(k)]] <- optRandomParC(M=facil_net,
                                                       k=k,
                                                       preSpecM=.50,
                                                       rep=ifelse(k==1,1,1000),
                                                       blocks=mcp(k),
                                                       approach="val",
                                                       nCores=ifelse(k==1,1,0))
}


#Saving these results so that we have something to work on later
saveRDS(res_trans_pos_bm,
        "results/res_trans_pos_bm.Rdata")


# Relative Fit Function
rf_trans_pos <- list()

for (k in 1:9){
  rf_trans_pos[[k]] <- RF(res = res_trans_pos_bm[[k]], 
                          m = 30, 
                          loops = TRUE)
}

rf_trans_pos

saveRDS(rf_trans_pos,
        "results/res_trans_pos_rf.Rdata")

# 3. Hierarchical -------------------------------------------------

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
  pos_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # turning everything to null
  pos_net_bm[1,,1:nbl] <- "nul"
  
  # 
  for (i in 1:(nbl-1)){
    pos_net_bm[1,(i+1),i] <- "com"
  }
  
  return(pos_net_bm)
}

printIMmodel<-function(IM)print(apply(IM,2:3,function(x)paste(na.omit(x),collapse=",")))
printIMmodel(mcp(10))

# fitting the blockmodel

res_hier_pos_bm <- list()

for (k in 2:10){
  res_hier_pos_bm[[as.character(k)]] <- optRandomParC(M=facil_net,
                                                      k=k,
                                                      preSpecM=.50,
                                                      rep=ifelse(k==1,1,1000),
                                                      blocks=mcp(k),
                                                      approach="val",
                                                      nCores=ifelse(k==1,1,0))
}

# Saving Results
saveRDS(res_hier_pos_bm,
        "results/res_hier_pos_bm.Rdata")


# Relative Fit Function
rf_hier_pos <- list()

for (k in 1:9){
  rf_hier_pos[[k]] <- RF(res = res_hier_pos_bm[[k]], 
                         m = 30, 
                         loops = TRUE)
}

rf_hier_pos

saveRDS(rf_hier_pos,
        "results/res_hier_pos_rf.Rdata")


# 4. Cohesive -------------------------------------------------

# Example with four blocks
#      1        2          3        4      
# 1   com      nul        nul      nul   
# 2   nul      com        nul      nul    
# 3   nul      nul        com      nul    
# 4   nul      nul        nul      com   


mcp <- function(nbl){
  pos_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # turning everything to null
  pos_net_bm[1,,1:nbl] <- "nul"
  
  for (i in 1:nbl){
    pos_net_bm[1,i,i] <- "com"
  }
  
  return(pos_net_bm)
}


# Fitting blockmodels
res_coh_pos_bm <- list()

for (k in 2:10){
  res_coh_pos_bm[[as.character(k)]] <- optRandomParC(M=facil_net,
                                                     k=k,
                                                     preSpecM=.50,
                                                     rep=ifelse(k==1,1,1000),
                                                     blocks=mcp(k),
                                                     approach="val",
                                                     nCores=ifelse(k==1,1,0))
}

#Saving these results 
saveRDS(res_coh_pos_bm,
        "results/res_coh_pos_bm.Rdata")


# Relative Fit Function
rf_coh_pos <- list()

for (k in 1:9){
  rf_coh_pos[[k]] <- RF(res = res_coh_pos_bm[[k]], 
                        m = 30, 
                        loops = TRUE)
}

rf_coh_pos

saveRDS(rf_coh_pos,
        "results/res_coh_pos_rf.Rdata")

# 5. Transitive Cohesive -------------------------------------------------

#same as transitive but now cohesive blocks
mcp <- function(nbl){
  pos_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # turning everything to null
  pos_net_bm[1,,1:nbl] <- "nul"
  
  for (i in 1:nbl){
    pos_net_bm[1,i,i] <- "com"
  }
  
  for (i in 1:(nbl-1)){
    pos_net_bm[1,(i+1):nbl,i] <- "com"
  }
  
  return(pos_net_bm)
}

printIMmodel(mcp(10))

# fitting the blockmodel

res_trans_coh_pos_bm <- list()

for (k in 2:10){
  res_trans_coh_pos_bm[[as.character(k)]] <- optRandomParC(M=facil_net,
                                                           k=k,
                                                           preSpecM=.50,
                                                           rep=ifelse(k==1,1,1000),
                                                           blocks=mcp(k),
                                                           approach="val",
                                                           nCores=ifelse(k==1,1,0))
}

res_trans_coh_pos_bm

#Saving these results 
saveRDS(res_trans_coh_pos_bm,
        "results/res_trans_coh_pos_bm.Rdata")


# Relative Fit Function
rf_trans_coh_pos <- list()

for (k in 1:9){
  rf_trans_coh_pos[[k]] <- RF(res = res_trans_coh_pos_bm[[k]], m = 30, loops = TRUE)
}

rf_trans_coh_pos

saveRDS(rf_trans_coh_pos,
        "results/res_trans_coh_pos_rf.Rdata")



# 6. Hierarchical Cohesve -------------------------------------------------

#same as hierarchical but now cohesive blocks

mcp <- function(nbl){
  pos_net_bm <- array(NA , dim = c(2, nbl, nbl))
  
  # turning everything to null
  pos_net_bm[1,,1:nbl] <- "nul"
  
  # 
  for (i in 1:nbl){
    pos_net_bm[1,i,i] <- "com"
  }
  
  for (j in 1:(nbl-1)){
    pos_net_bm[1,(j+1),j] <- "com"
  }
  
  return(pos_net_bm)
}

printIMmodel<-function(IM)print(apply(IM,2:3,function(x)paste(na.omit(x),collapse=",")))
printIMmodel(mcp(10))

# fitting the blockmodel
set.seed(2024)

res_hier_coh_pos_bm <- list()

for (k in 2:10){
  res_hier_coh_pos_bm[[as.character(k)]] <- optRandomParC(M=facil_net,
                                                          k=k,
                                                          preSpecM=.50,
                                                          rep=ifelse(k==1,1,1000),
                                                          blocks=mcp(k),
                                                          approach="val",
                                                          nCores=ifelse(k==1,1,0))
}

#Saving these results 
saveRDS(res_hier_coh_pos_bm,
        "results/res_hier_coh_pos_bm.Rdata")


# Relative Fit Function
rf_hier_coh_pos <- list()

for (k in 1:9){
  rf_hier_coh_pos[[k]] <- RF(res = res_hier_coh_pos_bm[[k]],
                             m = 30,
                             loops = TRUE)
}

rf_hier_coh_pos

saveRDS(rf_hier_coh_pos,
        "results/res_hier_coh_pos_rf.Rdata")


# The manuscript figures will be obtained in the next script


