library(tidyverse)
library(ggplot2)
library(cowplot)

df <- read.csv("/mnt/dv/wid/projects4/SolisLemus-network-merging/InPhyNet-Simulations/est-gts/data/out.csv") %>%
    mutate(
        retic_diff = as.factor(nretic_true - nretic_est),
        input_error = min_unrooted_nj_hwcd + sum_constraint_hwcd,
        ngt_label = paste0(ngt, " gt"),
        ils_label = paste0(str_to_title(ils), " ILS"),
        ntaxa_label = factor(
            paste0(ntaxa, " taxa"),
            levels = paste0(sort(unique(ntaxa)), " taxa")
        )
    )
nrow(df)

# Distance matrix estimation error:     min_unrooted_nj_hwcd
# Constraint network estimation error:  sum_constraint_hwcd
# Final inference error:                unrooted_hwcd



df %>%
    ggplot(aes(x = input_error, y = unrooted_hwcd, color = ils)) +
    facet_grid(ntaxa_label ~ ngt_label) +
    geom_jitter() +
    geom_abline(slope = 1, intercept = 0, color = "black", 
                linetype = "dashed", alpha = 0.5) +
    labs(x = "Input Error", y = "Output Error") +
    scale_color_manual("ILS", values = c('#1b9e77', '#d95f02', '#7570b3'))



diffs <- df$sum_constraint_hwcd + df$min_unrooted_nj_hwcd - df$unrooted_min_greedy_hwcd
qqnorm(diffs)
qqline(diffs)
