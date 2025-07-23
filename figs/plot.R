library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggdist)
library(gghalves)
library(ggh4x)

df <- read.csv("data/all.csv")
nrow(df)

theme_set(theme_bw())



df_plot <- df %>% 
  mutate(
    ntaxa_num = as.integer(as.character(ntaxa)),
    m = paste0("m = ", as.factor(m)),
    ntaxa_char = paste0(ntaxa, " tips"),
    imethod = if_else(imethod == "snaq", "SNaQ",
              if_else(imethod == "squirrel", "Squirrel",
              if_else(imethod == "phylonet", "PhyloNet-MPL",
              if_else(imethod == "phylonet-ml", "PhyloNet-ML", "NA"))))
  )
levels(df_plot$ntaxa_char) = c("30 tips", "50 tips", "100 tips", "200 tips")

# INPUT VS OUTPUT
p_inout <- ggplot(df_plot, aes(x = input_error, y = hwcd, color = ntaxa_char, shape = imethod)) +
    geom_jitter(width = 0.05, height = 0.1, stroke=0.5, size=0.85, alpha=0.6) +
    geom_abline(slope = 1, intercept = 0, color = "black", lty = "dashed") +
    scale_shape_manual(values = c("SNaQ" = 3, "Squirrel" = 1, "PhyloNet-MPL" = 8, "PhyloNet-ML" = 6)) +
    scale_color_manual(
      breaks = levels(df_plot$ntaxa_char),
      values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
    ) +
    labs(
      x = "Input Error (HWCD)",
      y = "Output Error (HWCD)",
      color = "Number of Tips",
      shape = "Method",
      size = ""
    ) +
    theme_bw() +
    scale_x_sqrt(limits = c(-0.05, 450)) +
    scale_y_sqrt(limits = c(-0.05, 450)) +
    theme(
      panel.grid.minor  = element_blank(),
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.spacing.y  = unit(-0.5, "lines")
    ) +
    guides(
      color = guide_legend(order = 2),
      shape = guide_legend(order = 1)
    )
p_inout

pdf("figs/accuracy/input-vs-output.pdf", width=5, height=5.5)
p_inout
dev.off()


df_plot <- df_plot %>%
  mutate(
    nbp_str = paste0(nbp, " base pairs"),
    ngt_str = paste0(ngt, " gene trees"),
    nbp = paste0("nbp = ", nbp),
    ngt = paste0("ngt = ", ngt)
  )



# <1e8 to fix some bugged outputs
p_rt <- filter(df_plot, runtime_serial < 1e8 & imethod == "SNaQ") %>% #filter(df_plot, runtime_serial < 1e8 & (imethod == "SNaQ" | imethod == "PhyloNet-ML" | imethod == "PhyloNet-MPL")) %>%
  mutate(runtime_serial = runtime_serial / 360) %>%
  ggplot(aes(x = ntaxa_num, y = runtime_serial)) +
  facet_grid(m ~ imethod, scales="free") +
  geom_violin(
    aes(group    = interaction(ntaxa_num, m)),
    position     = position_dodge(width = 15),
    width        = 15,
    color        = "black"
  ) +
  stat_smooth(
    aes(group = m),
    method   = "lm",
    formula  = y ~ x,
    se       = FALSE,
    size     = 0.5,
    linetype = "dashed",
    alpha    = 0.5,
    color    = "black"
  ) +
  scale_x_continuous(
    breaks = sort(unique(df_plot$ntaxa_num)),
    labels = sort(unique(df_plot$ntaxa))
  ) +
  labs(
    x = "Number of Tips",
    y = "Runtime (hours)"
  ) +
  theme_classic() +
  theme(
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = .4),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = .25),
    strip.text         = element_text(face = "bold", size = 8),
    legend.position    = "bottom"
  )
p_rt

pdf("figs/runtime/linear-runtime.pdf", width=5, height=5)
p_rt
dev.off()




##################
df_clean <- df %>% 
  mutate(
    ngt   = factor(ngt,  levels = sort(unique(ngt)),
                         labels  = paste0(sort(unique(ngt)))),
    nbp   = factor(nbp,  levels = sort(unique(nbp)),
                         labels  = paste0(sort(unique(nbp)))),
    m     = factor(m,    levels = sort(unique(m)),
                         labels  = paste0("m = ", sort(unique(m)))),
    ntaxa = factor(ntaxa,levels = sort(unique(ntaxa)),
                         labels  = paste0(sort(unique(ntaxa)), " tips")),
    imethod = factor(imethod,
                     levels = c("snaq","squirrel", "phylonet", "phylonet-ml"),
                     labels = c("SNaQ","SQUIRREL", "PhyloNet-MPL", "PhyloNet-ML")),
    ils = factor(ils, levels = c("high", "low"), labels = c("High", "Low"))
  )

p_hwcd <- ggplot(df_clean,
        aes(x = ils, y = hwcd, fill = nbp)) +
  geom_boxplot(width = .55, outlier.shape = NA, linewidth = .25) +
  #facet_grid(ntaxa ~ ngt + imethod + m, labeller = label_value, scales="free") +
  facet_nested(ntaxa ~ imethod + m, labeller = label_value, scales="free") +
  scale_fill_manual(
    values = c("100" = "#1b9e77", "1000" = "#d95f02"),
    name = "Number of Base Pairs"
  ) +
  labs(x = "Level of ILS",
       y = "Output Error (HWCD)") +
  theme_classic(base_size = 9) +
  theme(
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = .4),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = .25),
    strip.text         = element_text(face = "bold", size = 8),
    legend.position    = "bottom"
  ) +
  expand_limits(x = 0, y = 0)
p_hwcd


pdf("figs/accuracy/hwcd.pdf", width=5, height=5)
p_hwcd
dev.off()

