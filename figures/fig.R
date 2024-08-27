library(tidyverse)
library(ggplot2)

df <- read.csv("true-gts/output/n100.csv")

ggplot(df, aes(x = constraint_error_sum, y = HWCD, color = ils)) +
    geom_smooth(method = "lm", formula = y ~ x, se = F) +
    geom_jitter(alpha = 0.15, width = 1, height = 1) +
    facet_grid(ngt~max_subset_size)

