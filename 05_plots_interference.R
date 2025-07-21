# Set Up --------------------------------------------------------

# loading packages
library(blockmodeling)
library(sna)
library(tidyverse)

# Loading and Pre-Processing Blockmodeling data ------------------------------------------------------------

# Identity Interference Network 

interf_net <- readRDS("data/interf_net.RData")

# Co-occurrence matrix
n_matrix <- readRDS("data/n_matrix.RData")


# Abbreviating Names
rownames(interf_net)[rownames(interf_net)== "black or african american"] <- "Black"
rownames(interf_net)[rownames(interf_net)== "domestic violence survivor"] <- "domestic v. surv"
rownames(interf_net)[rownames(interf_net)== "member of a cultural or hobby activity organization"] <- "member cultural"
rownames(interf_net)[rownames(interf_net)== "member of a gym, fitness, sports or outdoor activity club"] <- "member fitness"
rownames(interf_net)[rownames(interf_net)== "person having a mental health diagnosis"] <- "p. mental health"
rownames(interf_net)[rownames(interf_net)== "person living with a disability"] <- "p. disability"
rownames(interf_net)[rownames(interf_net)== "person living with chronic illness/injury"] <- "p. chronic ill"
colnames(interf_net) <- rownames(interf_net)


# here making it so that identities that co-ccurred but 0%
# of respondents reported they facilitated each other
# appear as 0s and
# blanks appear as identities that never co-occurred in the matrices

interf_net2 <- interf_net
interf_net2[n_matrix==0] <- NA
interf_net2[interf_net2==0] <- 0.0001
interf_net2[n_matrix==0] <- 0

# Blockmodeling Results - Pre-specified model 

res_pres_neg_bm <- readRDS("results/res_pres_neg_bm.Rdata")

unlistcf <- function(x) unlist(sapply(x, function(y) min(y$err)))

neg_pre_cf_df <- data.frame(k = 2:10, 
                            cf = unlistcf(res_pres_neg_bm))

neg_pre_cf_df

# Density Matrix
res_pres_neg_dens <- readRDS("results/res_pres_neg_dens.Rdata") 

# Loading and Pre-processing Relative Fit Data ----------------------------

neg_rf_files <- c(
  neg_pre_rf       = "res_pres_neg_rf.Rdata",
  neg_trans_rf     = "res_trans_neg_rf.Rdata",
  neg_hier_rf      = "res_hier_neg_rf.Rdata",
  neg_coh_rf       = "res_coh_neg_rf.Rdata",
  neg_trans_coh_rf = "res_trans_coh_neg_rf.Rdata",
  neg_hier_coh_rf  = "res_hier_coh_neg_rf.Rdata"
)

# Common path
base_path <- "results/"

# Read files and extract rf values
neg_rf_list <- map(paste0(base_path, neg_rf_files), readRDS)

names(neg_rf_list) <- names(neg_rf_files)

extract_rf_df <- function(rf_object, name) {
  tibble(k = 2:10, !!name := unlist(map(rf_object, ~ .x$RF)))
}

rf_dfs <- imap(neg_rf_list, extract_rf_df)

# Combine all into one data frame
interf_rf_df <- reduce(rf_dfs, left_join, by = "k")

interf_rf_df <- interf_rf_df  |>  
  pivot_longer(ends_with("rf"), 
               values_to = "rf", 
               names_to = "method")

interf_rf_df <- interf_rf_df |>  
  mutate(Method = case_when(
    method == "neg_coh_rf" ~ "Cohesive", 
    method == "neg_hier_rf" ~ "Hierarchical", 
    method == "neg_hier_coh_rf" ~ "Hierarchical Cohesive", 
    method == "neg_pre_rf" ~ "Hypothesized", 
    method == "neg_trans_rf" ~ "Transitive", 
    method == "neg_trans_coh_rf" ~ "Transitive Cohesive")) %>%  
  mutate(Method = as.factor(Method),
         Method = fct_relevel(Method, 
                              c("Hypothesized", 
                                "Transitive", 
                                "Hierarchical", 
                                "Cohesive", 
                                "Transitive Cohesive", 
                                "Hierarchical Cohesive")))

interf_rf_df |>  
  filter(Method == "Hypothesized", 
         k == 5)

# Figure 5 - Relative Fit and Criterion Function --------------------------

# 1. Relative Fit - Figure 3 Panel A

png("figures/figure5_panela.png",
    width = 15,
    height = 10,
    units = "cm",
    res = 300)
ggplot(interf_rf_df, aes(x = k, y = rf, shape = Method, 
                      fill = Method)) +
  geom_line(color = "grey60")+
  geom_point(color = "grey20")+
  theme_minimal()+
  scale_shape_manual(values = c("Cohesive" = 15, "Hierarchical" = 17,
                                "Hierarchical Cohesive" = 8,
                                "Hypothesized" = 19,
                                "Transitive" = 13,
                                "Transitive Cohesive" = 3))+
  scale_y_continuous(limits = c(0,.5))+
  scale_x_continuous(breaks = c(0:10, 1))+
  labs(y = "RF") +
  geom_point(x = 5, y = 0.336, 
             shape = 19,
             size = 2,
             fill = "blue", 
             color = "blue")
dev.off()

# 2. Criterion Function - Figure 5 Panel B

png("figures/figure5_panelb.png",
    width = 10,
    height = 10,
    units = "cm",
    res = 300)
ggplot(neg_pre_cf_df, aes(x = k, y = cf)) +
  geom_line(color = "grey60")+
  geom_point(color = "grey20")+
  theme_minimal()+
  scale_y_continuous(limits = c(140,210))+
  scale_x_continuous(breaks = c(0:10, 1))+
  labs(y = "CF") +
  geom_point(x = 5, y = 173.1437,
             shape = 19,
             size = 2,
             fill = "blue", 
             color = "blue")
dev.off()


# Figure 6: Blockmodeling Results -----------------------------------------

neg_pre_dens3 <- funByBlocks(res_pres_neg_bm$'3')
neg_pre_dens5 <- funByBlocks(res_pres_neg_bm$'5')
neg_pre_dens7 <- funByBlocks(res_pres_neg_bm$'7')

neg_pre_dens3[is.nan(neg_pre_dens3)] <- 0
neg_pre_dens5[is.nan(neg_pre_dens5)] <- 0
neg_pre_dens7[is.nan(neg_pre_dens7)] <- 0

# Changing the order of the blocks
neg_pre_dens7 <-neg_pre_dens7[c(1,2,3,6,7,5,4),c(1,2,3,6,7,5,4)] 

rownames(neg_pre_dens7) <- 1:7
colnames(neg_pre_dens7) <- 1:7

new_clu <- data.frame(orig_clu = clu(res_pres_neg_bm[[5]])) |>  
  mutate(new_clu = case_when(
    orig_clu == 1 ~ 1,
    orig_clu == 2 ~ 2,
    orig_clu == 3 ~ 3,
    orig_clu == 4 ~ 7,
    orig_clu == 5 ~ 6,
    orig_clu == 6 ~ 4, 
    orig_clu == 7 ~ 5))

png("figures/figure6.png",
    width = 13,
    height = 15,
    units = "cm",
    res = 300)
par(mfrow = c(3, 2))
plotMat(interf_net2, 
        clu = clu(res_pres_neg_bm$'3'), 
        print.digits.cells = 3, 
        main = "",
        mar = c(1,1,3,1), 
        printMultipliedMessage = F, 
        #print.val = F, 
        plot.legend = F)
plotMat(neg_pre_dens3, 
        #print.digits.cells = 3,
        main = "", 
        mar = c(1,1,1,1), 
        printMultipliedMessage = F, 
        print.val = T, 
        plot.legend = F)
plotMat(interf_net2, 
        clu = clu(res_pres_neg_bm$'5'), 
        print.digits.cells = 3, 
        main = "", 
        mar = c(1,1,2,1), 
        printMultipliedMessage = F, 
        #print.val = F, 
        plot.legend = F)
plotMat(neg_pre_dens5, 
        #print.digits.cells = 3,
        main = "", 
        mar = c(1,1,1,1), 
        printMultipliedMessage = F, 
        print.val = T, 
        plot.legend = F)
plotMat(interf_net2, 
        clu = new_clu$new_clu, 
        print.digits.cells = 3, 
        main = "", 
        mar = c(1,1,2,1), 
        printMultipliedMessage = F, 
        #print.val = F, 
        plot.legend = F)
plotMat(neg_pre_dens7, 
        #print.digits.cells = 3,
        main = "", 
        mar = c(1,1,1,1), 
        printMultipliedMessage = F, 
        print.val = T, 
        plot.legend = F)
dev.off()

