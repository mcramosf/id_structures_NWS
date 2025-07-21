
# Load packages --------------------------------------------------------
library(blockmodeling)
library(sna)
library(tidyverse)

# Data Import -------------------------------------------------------------

# Co-occurrence Matrix

# The network contains 37 identities. 
# Each cell indicates the number of respondents
# that simultaneously held the corresponding pair of identities. 

n_matrix <- readRDS("data/n_matrix.RData")

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

# Figure 2 - Co-occurrence Matrix ----------------------------------------------------------------

# Abbreviating names
n_matrix_plot <- n_matrix

rownames(n_matrix_plot)[rownames(n_matrix_plot)== "black or african american"] <- "Black"
rownames(n_matrix_plot)[rownames(n_matrix_plot)== "domestic violence survivor"] <- "domestic v. surv"
rownames(n_matrix_plot)[rownames(n_matrix_plot)== "member of a cultural or hobby activity organization"] <- "member cultural"
rownames(n_matrix_plot)[rownames(n_matrix_plot)== "member of a gym, fitness, sports or outdoor activity club"] <- "member fitness"
rownames(n_matrix_plot)[rownames(n_matrix_plot)== "person having a mental health diagnosis"] <- "p. mental health"
rownames(n_matrix_plot)[rownames(n_matrix_plot)== "person living with a disability"] <- "p. disability"
rownames(n_matrix_plot)[rownames(n_matrix_plot)== "person living with chronic illness/injury"] <- "p. chronic ill"
colnames(n_matrix_plot) <- rownames(n_matrix_plot)

# Producing the plot
png("figures/figure2_co_occurrence.png",
    width = 10,
    height = 10,
    units = "cm",
    res = 300)
plotMat(n_matrix_plot, 
        print.val = TRUE, 
        print.0 = T,
        print.digits.cells = 3, 
        main = "")
dev.off()

# Descriptives presented for the co-ocurrence matrix

# obtaining the number of values = 0 and median cell value

n_df <- as.data.frame(as.table(n_matrix))

n_df <- n_df %>% 
  mutate(Var1 = as.character(Var1), 
         Var2 = as.character(Var2)) %>% 
  filter(Var1 != Var2) %>% 
  mutate(newvar = NA) 

for (i in 1:nrow(n_df)){
  
  n_df[i,"newvar"] <-str_c(sort(c(n_df[i,"Var1"], 
                                  n_df[i,"Var2"]))[1], 
                           sort(c(n_df[i,"Var1"], 
                                  n_df[i,"Var2"]))[2])
}

n_df <- n_df %>% 
  distinct(across(newvar), 
           .keep_all = T)
head(n_df)

# 61 cell values = 0
table(n_df$Freq)

# median cell value = 8 
summary(n_df$Freq)

# exploring in more detail which cells have a cell
# value  = 0
n_df0 <- n_df %>% 
  filter(Freq == 0)

male <- c("son", 
          "male", 
          "husband", 
          "father", 
          "brother",
          "boyfriend")

female <- c("wife", 
            "sister", 
            "mother", 
            "girlfriend", 
            "female", 
            "daughter")

n_df_gendered <- n_df %>% 
  filter((Var1 %in% male |
            Var2 %in% male) &
           (Var1 %in% female |
              Var2 %in% female))

table(n_df_gendered$Freq)

# n df without gendered roles

n_df_non_gendered <- n_df %>% 
  anti_join(n_df_gendered)

filter(n_df_non_gendered, 
       n_df_non_gendered$Freq == 0)

# 11 of these are mutually exclusive identities
# others are very unlikely to occur (atheist and churchgoer)


# Descriptive Statistics Facilitation Networks ----------------------------

# Density 
# number of ties greater than 0 out of possible ties
# excludes the diagonal
facil_df <- as.data.frame(as.table(facil_net))
sum(facil_net>0)
1011/(37^2-37)

facil_df %>% 
  filter(Freq >0) %>% 
  count()
1011/(666*2)

# Average tie value given that tie is non-zero
summary(facil_df$Freq[facil_df$Freq>0])


# Descriptive Statistics Interfering Networks  ------------------------------------------------

interf_df <- as.data.frame(as.table(interf_net))

# Density 
# number of ties greater than 0 out of possible ties
# excludes the diagonal
interf_df %>% 
  filter(Freq >0) %>% 
  count()
884/(666*2)

# Average tie value given that tie is non-zero
summary(interf_df$Freq[interf_df$Freq>0])
