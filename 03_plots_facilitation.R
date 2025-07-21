# Set Up --------------------------------------------------------

# loading packages
library(blockmodeling)
library(sna)
library(tidyverse)

# Loading and Pre-Processing Blockmodeling data ------------------------------------------------------------

# Identity Facilitation Network 

facil_net <- readRDS("data/facil_net.RData")

# Co-occurrence matrix
n_matrix <- readRDS("data/n_matrix.RData")

# here making it so that identities that co-ccurred but 0%
# of respondents reported they facilitated each other
# appear as 0s and
# blanks appear as identities that never co-occurred in the matrices

facil_net2 <- facil_net
facil_net2[n_matrix==0] <- NA
facil_net2[facil_net2==0] <- 0.0001
facil_net2[n_matrix==0] <- 0

# Abbreviating Names
rownames(facil_net)[rownames(facil_net)== "black or african american"] <- "Black"
rownames(facil_net)[rownames(facil_net)== "domestic violence survivor"] <- "domestic v. surv"
rownames(facil_net)[rownames(facil_net)== "member of a cultural or hobby activity organization"] <- "member cultural"
rownames(facil_net)[rownames(facil_net)== "member of a gym, fitness, sports or outdoor activity club"] <- "member fitness"
rownames(facil_net)[rownames(facil_net)== "person having a mental health diagnosis"] <- "p. mental health"
rownames(facil_net)[rownames(facil_net)== "person living with a disability"] <- "p. disability"
rownames(facil_net)[rownames(facil_net)== "person living with chronic illness/injury"] <- "p. chronic ill"
colnames(facil_net) <- rownames(facil_net)

# Blockmodeling Results - Pre-specified model 

res_pres_pos_bm <- readRDS("results/res_pres_pos_bm.Rdata")

unlistcf <- function(x) unlist(sapply(x, function(y) min(y$err)))

pos_pre_cf_df <- data.frame(k = 2:10, 
                            cf = unlistcf(res_pres_pos_bm))

pos_pre_cf_df

# Density Matrix
res_pres_pos_dens <- readRDS("results/res_pres_pos_dens.Rdata") 

res_pres_pos_dens[[5]]

# Loading and Pre-processing Relative Fit Data ----------------------------

pos_rf_files <- c(
  pos_pre_rf       = "res_pres_pos_rf.Rdata",
  pos_trans_rf     = "res_trans_pos_rf.Rdata",
  pos_hier_rf      = "res_hier_pos_rf.Rdata",
  pos_coh_rf       = "res_coh_pos_rf.Rdata",
  pos_trans_coh_rf = "res_trans_coh_pos_rf.Rdata",
  pos_hier_coh_rf  = "res_hier_coh_pos_rf.Rdata"
)

# Common path
base_path <- "results/"

# Read files and extract rf values
pos_rf_list <- map(paste0(base_path, pos_rf_files), readRDS)

names(pos_rf_list) <- names(pos_rf_files)

extract_rf_df <- function(rf_object, name) {
  tibble(k = 2:10, !!name := unlist(map(rf_object, ~ .x$RF)))
}

rf_dfs <- imap(pos_rf_list, extract_rf_df)

# Combine all into one data frame
facil_rf_df <- reduce(rf_dfs, left_join, by = "k")

facil_rf_df <- facil_rf_df  |>  
  pivot_longer(ends_with("rf"), 
               values_to = "rf", 
               names_to = "method")

facil_rf_df <- facil_rf_df |>  
  mutate(Method = case_when(
    method == "pos_coh_rf" ~ "Cohesive", 
    method == "pos_hier_rf" ~ "Hierarchical", 
    method == "pos_hier_coh_rf" ~ "Hierarchical Cohesive", 
    method == "pos_pre_rf" ~ "Hypothesized", 
    method == "pos_trans_rf" ~ "Transitive", 
    method == "pos_trans_coh_rf" ~ "Transitive Cohesive")) %>%  
  mutate(Method = as.factor(Method),
         Method = fct_relevel(Method, 
                              c("Hypothesized", 
                                "Transitive", 
                                "Hierarchical", 
                                "Cohesive", 
                                "Transitive Cohesive", 
                                "Hierarchical Cohesive")))

facil_rf_df |>  
  filter(Method == "Hypothesized", 
         k == 4)

# Figure 3 - Relative Fit and Criterion Function --------------------------

# 1. Relative Fit - Figure 3 Panel A

png("figures/figure3_panela.png",
    width = 15,
    height = 10,
    units = "cm",
    res = 300)
ggplot(facil_rf_df, aes(x = k, y = rf, shape = Method, 
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
  geom_point(x = 4, y = 0.246, 
             shape = 19,
             size = 2,
             fill = "blue", 
             color = "blue")
dev.off()

# 2. Criterion Function - Figure 3 Panel B

png("figures/figure3_panelb.png",
    width = 10,
    height = 10,
    units = "cm",
    res = 300)
ggplot(pos_pre_cf_df, aes(x = k, y = cf)) +
  geom_line(color = "grey60")+
  geom_point(color = "grey20")+
  theme_minimal()+
  scale_y_continuous(limits = c(160,260))+
  scale_x_continuous(breaks = c(0:10, 1))+
  labs(y = "CF") +
  geom_point(x = 4, y = 211.1040,
             shape = 19,
             size = 2,
             fill = "blue", 
             color = "blue")
dev.off()

# 3. Blockmodeling Results - Figure 4

# Changing the order of the blocks

res_pres_pos_dens2 <- res_pres_pos_dens[[5]]
res_pres_pos_dens2 <- res_pres_pos_dens2[c(1,2,4,5,6,3),c(1,2,4,5,6,3)]

new_clu <- data.frame(orig_clu = clu(res_pres_pos_bm[[5]])) |>  
  mutate(new_clu = case_when(
    orig_clu == 1 ~ 1,
    orig_clu == 2 ~ 4,
    orig_clu == 3 ~ 2,
    orig_clu == 4 ~ 6,
    orig_clu == 5 ~ 3,
    orig_clu == 6 ~ 5))

png("figures/figure4.png",
    width = 25,
    height = 12,
    units = "cm",
    res = 300)
par(mfrow = c(1, 2))
plotMat(facil_net2, 
        clu = new_clu$new_clu, 
        print.digits.cells = 3, 
        main = "", 
        mar = c(2,1,3,1), 
        printMultipliedMessage = F, 
        #print.val = F, 
        plot.legend = F)
plotMat(res_pres_pos_dens2, 
        #print.digits.cells = 3,
        main = "", 
        mar = c(1,1,1,1), 
        printMultipliedMessage = F, 
        print.val = T, 
        plot.legend = F)
dev.off()
