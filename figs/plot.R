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
    ntaxa_char = paste0(ntaxa, " taxa"),
    imethod = if_else(imethod == "snaq", "SNaQ",
              if_else(imethod == "squirrel", "Squirrel",
              if_else(imethod == "phylonet", "PhyloNet-MPL",
              if_else(imethod == "phylonet-ml", "PhyloNet-ML", "NA"))))
  )
levels(df_plot$ntaxa_char) = c("25 taxa", "50 taxa", "100 taxa", "200 taxa")

# INPUT VS OUTPUT
p_inout <- df_plot %>%
  mutate(
    hwcd = hwcd + if_else(hwcd == 0, 0, runif(n(), -1, 1)),
    input_error = input_error + if_else(input_error == 0, 0, runif(n(), -1, 1))
  ) %>%
  ggplot(aes(x = input_error, y = hwcd, color = ntaxa_char, shape = imethod)) +
    #geom_jitter(width = 0.05, height = 0.1, stroke=0.5, size=0.85, alpha=0.6) +
    geom_point(size = 0.55, alpha = 0.6, stroke = 0.3) +
    geom_abline(slope = 1, intercept = 0, color = "black", lty = "dashed") +
    scale_shape_manual(values = c("SNaQ" = 3, "Squirrel" = 1, "PhyloNet-MPL" = 8, "PhyloNet-ML" = 6)) +
    scale_color_manual(
      breaks = levels(df_plot$ntaxa_char),
      values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
    ) +
    labs(
      x = "Input Error (HWCD)",
      y = "Output Error (HWCD)",
      color = "Number of Taxa",
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
    ngt = paste0("ngt = ", ngt),
    runtime_serial = runtime_serial / 360,
    ils = factor(ils, levels = c("low", "high"), labels = c("Low ILS", "High ILS"))
  )



# <1e8 to fix some bugged outputs
#p_rt <- filter(df_plot, runtime_serial < 1e8 & imethod == "SNaQ") %>%
p_rt <- df_plot %>%#filter(df_plot, runtime_serial < 1e8 & m == "m = 10" & imethod == "SNaQ") %>%
  ggplot(aes(x = ntaxa_num, y = runtime_serial, fill = as.factor(ils))) +
  facet_grid(imethod ~ m, scales="free") +
  geom_boxplot(
    aes(group    = interaction(ntaxa_num, m, ils)),
    position     = position_dodge(width = 5),
    width        = 5,
    color        = "black"
  ) +
  stat_smooth(
    aes(group = interaction(ntaxa_num, m, ils)),
    method   = "lm",
    formula  = y ~ x,
    se       = FALSE,
    linewidth= 0.5,
    linetype = "dashed",
    alpha    = 0.5,
    color    = "black"
  ) +
  scale_x_continuous(
    breaks = sort(unique(df_plot$ntaxa_num)),
    labels = sort(unique(df_plot$ntaxa))
  ) +
  labs(
    x = "Number of Taxa",
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




# rt_plot <- function(dat, method) {
#   dat <- dat %>%
#     filter(imethod == method & ils == "Low ILS")
#   summary_df <- dat %>%
#     group_by(m, ils, ntaxa_num) %>%
#     summarise(
#       mean_rt = median(runtime_serial),
#       ymin = quantile(runtime_serial, 0.1),
#       ymax = quantile(runtime_serial, 0.9)
#     )
#   dat %>%
#     ggplot(aes(x = ntaxa_num, y = runtime_serial)) +
#     facet_grid(m ~ ils, scales = "free") +
#     geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", linetype = "dashed") +
#     geom_point(data=summary_df, aes(x = ntaxa_num, y = mean_rt)) +
#     geom_errorbar(data=summary_df, aes(ymin = ymin, ymax = ymax, y = mean_rt)) +
#     labs(
#       x = "Number of Taxa",
#       y = "Runtime (hours)"
#     ) +
#     theme_classic() +
#     theme(
#       panel.border       = element_rect(colour = "black", fill = NA, linewidth = .4),
#       panel.grid.major.x = element_blank(),
#       panel.grid.major.y = element_line(colour = "grey90", linewidth = .25),
#       strip.text         = element_text(face = "bold", size = 8),
#       legend.position    = "bottom"
#     ) +
#     scale_x_continuous(
#       breaks = unique(df$ntaxa),
#       labels = unique(df$ntaxa)
#     ) +
#     ggtitle(method)
# }
# rt_plot(df_plot, "PhyloNet-MPL") + rt_plot(df_plot, "SNaQ") + rt_plot(df_plot, "Squirrel")

df_rt <- df_plot %>%
  filter(ils == "Low ILS")

summary_df <- df_rt %>%
  filter(imethod != "Squirrel" & imethod != "PhyloNet-ML") %>%
  group_by(m, ils, ntaxa_num, imethod) %>%
  summarise(
    mean_rt = mean(runtime_serial),
    ymin = quantile(runtime_serial, 0.025),
    ymax = quantile(runtime_serial, 0.975)
  )
nosquirrelnophyloml <- df_rt %>%
  filter(imethod != "Squirrel" & imethod != "PhyloNet-ML") %>%
  ggplot(aes(x = ntaxa_num, y = runtime_serial)) +
  facet_grid(m ~ imethod, scales = "free") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", linetype = "dashed") +
  geom_point(data=summary_df, aes(x = ntaxa_num, y = mean_rt)) +
  geom_errorbar(data=summary_df, aes(ymin = ymin, ymax = ymax, y = mean_rt)) +
  labs(
    x = "Number of Taxa",
    y = "Runtime (hours)"
  ) +
  theme_classic() +
  theme(
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = .4),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = .25),
    strip.text         = element_text(face = "bold", size = 8),
    legend.position    = "bottom"
  ) +
  scale_x_continuous(
    breaks = unique(df$ntaxa),
    labels = unique(df$ntaxa)
  )
summary_df <- df_rt %>%
  filter(imethod == "Squirrel") %>%
  group_by(m, ils, ntaxa_num, imethod) %>%
  summarise(
    mean_rt = mean(runtime_serial),
    ymin = quantile(runtime_serial, 0.025),
    ymax = quantile(runtime_serial, 0.975)
  )
squirrel <- df_rt %>%
  filter(imethod == "Squirrel") %>%
  ggplot(aes(x = ntaxa_num, y = runtime_serial)) +
  facet_grid(m ~ imethod, scales = "free") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", linetype = "dashed") +
  geom_point(data=summary_df, aes(x = ntaxa_num, y = mean_rt)) +
  geom_errorbar(data=summary_df, aes(ymin = ymin, ymax = ymax, y = mean_rt)) +
  labs(
    x = "Number of Taxa",
    y = "Runtime (hours)"
  ) +
  theme_classic() +
  theme(
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = .4),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = .25),
    strip.text         = element_text(face = "bold", size = 8),
    legend.position    = "bottom"
  ) +
  scale_x_continuous(
    breaks = unique(df$ntaxa),
    labels = unique(df$ntaxa)
  ) +
  ylim(0, NA)
summary_df <- df_rt %>%
  filter(imethod == "PhyloNet-ML") %>%
  group_by(m, ils, ntaxa_num, imethod) %>%
  summarise(
    mean_rt = mean(runtime_serial),
    ymin = quantile(runtime_serial, 0.025),
    ymax = quantile(runtime_serial, 0.975)
  )
phylonetml <- df_rt %>%
  filter(imethod == "PhyloNet-ML") %>%
  ggplot(aes(x = ntaxa_num, y = runtime_serial)) +
  facet_grid(m ~ imethod, scales = "free") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red", linetype = "dashed") +
  geom_point(data=summary_df, aes(x = ntaxa_num, y = mean_rt)) +
  geom_errorbar(data=summary_df, aes(ymin = ymin, ymax = ymax, y = mean_rt)) +
  labs(
    x = "Number of Taxa",
    y = ""
  ) +
  theme_classic() +
  theme(
    panel.border       = element_rect(colour = "black", fill = NA, linewidth = .4),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = .25),
    strip.text         = element_text(face = "bold", size = 8),
    legend.position    = "bottom"
  ) +
  scale_x_continuous(
    breaks = unique(df$ntaxa),
    labels = unique(df$ntaxa)
  ) +
  ylim(0, NA)

p_rt <- (squirrel + xlab("")) + (nosquirrelnophyloml + ylab("")) + ((phylonetml + xlab("") + ylab("")) / plot_spacer()) +
  plot_layout(widths = c(1, 2.05, 1), axis_titles = "collect")
p_rt


pdf("figs/runtime/linear-methods-separated.pdf", width=11, height=5)
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
                         labels  = paste0(sort(unique(ntaxa)), " taxa")),
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






gen_plot <- function() {
  df_plot %>%
    mutate(ninput = factor(nbp * ngt), nbp = factor(nbp), ngt = factor(ngt)) %>%
    mutate(ngt = paste0("ngt = ", ngt), nbp = paste0("nbp = ", nbp)) %>%
    mutate(
      hwcd = hwcd + if_else(hwcd == 0, 0, runif(n(), -1, 1)),
      input_error = input_error + if_else(input_error == 0, 0, runif(n(), -1, 1)),
      ils = factor(paste0(str_to_title(ils), " ILS"), levels = c("Low ILS", "High ILS"))
    ) %>%
    ggplot(aes(x = input_error, y = hwcd, color = ntaxa_char, shape = imethod)) +
      geom_point(size = 0.35, alpha = 0.45, stroke = 0.3) +
      geom_abline(slope = 1, intercept = 0, color = "black", lty = "dashed", alpha = 0.8) +
      scale_shape_manual(values = c("SNaQ" = 3, "Squirrel" = 1, "PhyloNet-MPL" = 8, "PhyloNet-ML" = 6)) +
      scale_color_manual(
        breaks = levels(df_plot$ntaxa_char),
        values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
      ) +
      labs(
        x = "Input Error (HWCD)",
        y = "Output Error (HWCD)",
        color = "Number of Taxa",
        shape = "Method",
        size = ""
      ) +
      theme_bw() +
      scale_x_sqrt(limits = c(0.0, 450)) +
      scale_y_sqrt(limits = c(0.0, 450)) +
      expand_limits(x = 0, y = 0) +
      theme_classic(base_size = 9) +
      theme(
        panel.border       = element_rect(colour = "black", fill = NA, linewidth = .4),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey90", linewidth = .25),
        strip.text         = element_text(face = "bold", size = 8),
        legend.position    = "bottom",
        legend.box         = "vertical",
        legend.spacing.y   = unit(-0.5, "lines")
      ) +
      guides(
        color = guide_legend(
          order = 2,
          override.aes = list(size = 1.5, alpha = 1.0)
        ),
        shape = guide_legend(
          order = 1,
          override.aes = list(size = 1.5, alpha = 1.0)
        )
      ) +
      facet_grid(ils ~ m)
      #facet_nested(ils + m ~ ngt + nbp)
}
df_plot <- df %>% 
  mutate(
    ntaxa_num = as.integer(as.character(ntaxa)),
    m = paste0("m = ", as.factor(m)),
    ntaxa_char = paste0(ntaxa, " taxa"),
    imethod = if_else(imethod == "snaq", "SNaQ",
              if_else(imethod == "squirrel", "Squirrel",
              if_else(imethod == "phylonet", "PhyloNet-MPL",
              if_else(imethod == "phylonet-ml", "PhyloNet-ML", "NA"))))
  )
levels(df_plot$ntaxa_char) = c("25 taxa", "50 taxa", "100 taxa", "200 taxa")

pdf("figs/accuracy/input-vs-output.pdf", width=5, height=5.5)
gen_plot()
dev.off()


ggplot(df_plot, aes(x = ngt, y = hwcd - input_error)) +
  geom_point() +
  facet_grid(ntaxa ~ nbp)
