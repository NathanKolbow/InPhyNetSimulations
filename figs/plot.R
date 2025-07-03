library(ggplot2)
library(tidyverse)

df <- read.csv("data/all.csv")

ggplot(df, aes(x = input_error, y = hwcd, color = imethod, shape = as.factor(ngt), size = as.factor(nbp))) +
    facet_grid(m~ntaxa) +
    geom_jitter(width = 0.15, height = 0) +
    geom_abline(slope = 1, intercept = 0, color = "black", lty = "dashed") +
    scale_size_manual(values = c("100" = 1, "1000" = 3)) +
    scale_shape_manual(values = c("100" = 3, "1000" = 1))


ggplot(df, aes(x = as.factor(ngt), y = hwcd, color = as.factor(nbp))) +
    facet_grid(ntaxa~imethod) +
    geom_boxplot()
