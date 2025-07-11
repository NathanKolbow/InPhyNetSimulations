library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggdist)
library(gghalves)

df <- read.csv("data/all.csv")
nrow(df)

ggplot(df, aes(x = input_error, y = hwcd, color = imethod, shape = as.factor(ngt), size = as.factor(nbp))) +
    facet_grid(m~ntaxa) +
    geom_jitter(width = 0.05, height = 0.05) +
    geom_abline(slope = 1, intercept = 0, color = "black", lty = "dashed") +
    scale_size_manual(values = c("100" = 1, "1000" = 3)) +
    scale_shape_manual(values = c("100" = 3, "1000" = 1))


ggplot(df, aes(x = as.factor(ntaxa), y = runtime_serial, fill = as.factor(ngt))) +
  facet_grid(imethod ~ nbp, scale="free") +
  geom_boxplot()







##################
df_clean <- df %>% 
  mutate(
    ngt   = factor(ngt,  levels = sort(unique(ngt)),
                         labels  = paste0("ngt = ", sort(unique(ngt)))),
    nbp   = factor(nbp,  levels = sort(unique(nbp)),
                         labels  = paste0("nbp = ", sort(unique(nbp)))),
    m     = factor(m,    levels = sort(unique(m)),
                         labels  = paste0("m = ", sort(unique(m)))),
    ntaxa = factor(ntaxa,levels = sort(unique(ntaxa)),
                         labels  = paste0("ntaxa = ", sort(unique(ntaxa)))),
    imethod = factor(imethod,
                     levels = c("snaq","squirrel"),
                     labels = c("SNaQ","SQUIRREL"))
  )
df_clean$hwcd <- df_clean$runtime_serial
p_hwcd <- ggplot(df_clean,
        aes(x = ngt, y = hwcd, fill = nbp)) +

  geom_boxplot(width = .55, outlier.shape = NA, linewidth = .25) +

  facet_grid(ntaxa ~ imethod + m, labeller = label_value, scales="free") +

  scale_fill_brewer(palette = "Dark2", name = "Number of base pairs (nbp)") +

  labs(x = "Number of gene trees (ngt)",
       y = "HWCD") +

  theme_classic(base_size = 9) +
  theme(
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = .4),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = .25),
    strip.text         = element_text(face = "bold", size = 8),
    legend.position    = "top"
  )
p_hwcd
